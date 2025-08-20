
import argparse
import numpy as np
import pandas as pd
import torch
import pyro
import pyro.distributions as dist
from pyro.infer import SVI, Trace_ELBO, infer_discrete
from pyro.infer.autoguide import AutoDelta
from pyro.optim import Adam
from pyro import poutine
from sklearn.ensemble import IsolationForest
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde




def parse_args():
    parser = argparse.ArgumentParser(description="Fetal Fraction Estimation and Genotyping from cfDNA")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file (cfDNA variants)")
    parser.add_argument("-o", "--output", default="genotype_output.csv", help="Output CSV file for genotype probabilities")
    parser.add_argument("--nsteps", type=int, default=2000, help="Number of SVI steps")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility")
    return parser.parse_args()


def load_and_filter_data(file_path, seed):
    df = pd.read_csv(file_path , sep="\t")
    vaf_col='vaf'
    frag_col='frag_stat'
    df = df[(df[vaf_col].between(0.025, 0.975)) & (df[frag_col].between(-4,4))]
    #clf = IsolationForest(contamination=0.05, random_state=seed)
    #df = df[clf.fit_predict(df[['VAF', 'frag_stat']]) == 1]
    #return df
    df['vaf'] = df['vaf'].replace(r'^\s*$', np.nan, regex=True) 
    df['frag_stat'] = df['frag_stat'].replace(r'^\s*$', np.nan, regex=True)    # Remplace les espaces par NaN
    df['vaf'] = df['vaf'].fillna(df['vaf'].median())    
    df['frag_stat'] = df['frag_stat'].fillna(df['frag_stat'].median())    
   # Remplace NaN par la médiane
    for col in [vaf_col,frag_col]:
      df[col] = pd.to_numeric(df[col], errors='coerce')
      
    return df


def estimate_initial_ff(df):
    vaf_vals = df['vaf'].values
    kde = gaussian_kde(vaf_vals, bw_method=0.2)
    x_eval = np.linspace(0.025, 0.975, 1000)
    density = kde(x_eval)
    local_max_idx = argrelextrema(density, np.greater)[0]
    peak_vafs = x_eval[local_max_idx]
    if len(local_max_idx) == 0:
      print("No local maxima found in VAF density — using fallback median-based estimate.")
      f_init = 2 * np.median(vaf_vals[vaf_vals > 0.5]) if np.any(vaf_vals > 0.5) else 0.1
    else:
      peak_vafs = x_eval[local_max_idx]
      f_init = 2 * (1 - peak_vafs[np.argmax(density[local_max_idx])])
      print(f"Estimated initial fetal fraction: {f_init:.3f}")
    return f_init


def initialize_frag_means(df, f_init):
    df['dist_to_cluster4_vaf'] = np.abs(df['vaf'] - (1 - f_init / 2))
    cluster4_sites = df.nsmallest(500, 'dist_to_cluster4_vaf')
    frag_stat_mean = cluster4_sites['frag_stat'].median()
    return frag_stat_mean * np.array([-1.0, 0.5, 0.0, -0.5, 1.0])


def build_model(X, frag_means):
    N, K = X.shape[0], 5

    def model(X):
        f = pyro.sample("f", dist.Uniform(0.01, 0.5))

        vaf_means = torch.stack([
            f / 2, (1 - f) / 2, torch.tensor(0.5), f + (1 - f) / 2, 1 - f / 2
        ])  # shape: [K]

        frag_means_ = pyro.param("frag_means", frag_means, constraint=dist.constraints.real)  # shape: [K]

        # Combine into 2D cluster means: shape [K, 2]
        cluster_means = torch.stack([vaf_means, frag_means_], dim=1)  # shape: [K, 2]

        # Sample standard deviations per cluster, per dimension (vaf and frag): shape [K, 2]
        vaf_stds = pyro.sample("vaf_stds", dist.LogNormal(torch.zeros(K), 0.2).to_event(1))
        frag_stds = pyro.sample("frag_stds", dist.LogNormal(torch.zeros(K), 0.2).to_event(1))  # [K]
        cluster_stds = torch.stack([vaf_stds, frag_stds], dim=1)  # shape [K, 2]

        # Construct diagonal covariance matrices: shape [K, 2, 2]
        cov_matrices = torch.diag_embed(cluster_stds ** 2)

        # Mixture weights
        weights = pyro.sample("weights", dist.Dirichlet(torch.ones(K)))
         
        with pyro.plate("data", N):
            assignment = pyro.sample("assignment",dist.Categorical(weights),infer={"enumerate": "parallel"})
            pyro.sample(
                "obs",
                dist.MultivariateNormal(cluster_means[assignment], cov_matrices[assignment]),
                obs=X
            )

    return model




def run_inference(X, frag_means, n_steps):
    model = build_model(X, frag_means)

    # Block discrete sites for autoguide and enumerate assignment
    guide = AutoDelta(poutine.block(model, hide=["assignment"]))
    
    svi = SVI(model, guide, Adam({"lr": 0.02}), loss=Trace_ELBO())

    pyro.clear_param_store()
    for step in range(n_steps):
        loss = svi.step(X)
        if step % 100 == 0:
            print(f"[{step}] Loss: {loss:.2f}")
    print(guide)
    return guide


def compute_posterior(X, guide, frag_means):
    if guide is None:
        raise ValueError("guide is None — make sure you've run inference and passed the trained guide.")

    N = X.shape[0]

    # Extract MAP values from the guide
    guide_trace = guide()
    f_map = guide_trace["f"]
    print(f"the fetal fraction is :{f_map}")
    weights_map = guide_trace["weights"]
    vaf_stds_map = guide_trace["vaf_stds"]
    frag_stds_map = guide_trace["frag_stds"]

    # Get frag means (learned parameter)
    frag_means_map = pyro.param("frag_means").detach()

    # Compute vaf means based on MAP f
    vaf_means_map = torch.stack([
        f_map / 2,
        (1 - f_map) / 2,
        torch.tensor(0.5, device=f_map.device),
        f_map + (1 - f_map) / 2,
        1 - f_map / 2])
        
    
    # Now compute posterior responsibilities
    # Build mixture components: shape [K, 2]
    means = torch.stack([vaf_means_map, frag_means_map], dim=1)
    stds = torch.stack([vaf_stds_map, frag_stds_map], dim=1)
    covs = torch.diag_embed(stds ** 2)

    # Compute per-cluster likelihoods
    component_dists = dist.MultivariateNormal(means, covs)
    log_probs = component_dists.log_prob(X[:, None, :])  # [N, K]

    # Add log weights
    log_weighted_probs = log_probs + torch.log(weights_map)

    # Normalize to get posterior probabilities
    posterior_log_probs = log_weighted_probs - torch.logsumexp(log_weighted_probs, dim=1, keepdim=True)
    posterior_probs = torch.exp(posterior_log_probs)  # shape [N, K]
    print(posterior_probs)
    return posterior_probs


def compute_genotype_probs(df, posterior_probs):
     
    fetal_00 = posterior_probs[:, 1]
    fetal_01 = posterior_probs[:, 0] + posterior_probs[:, 2] + posterior_probs[:, 4]
    fetal_11 = posterior_probs[:, 3]

    maternal_00 = posterior_probs[:, 0]
    maternal_01 = posterior_probs[:, 1] + posterior_probs[:, 2] + posterior_probs[:, 3]
    maternal_11 = posterior_probs[:, 4]

    genotype_df = pd.DataFrame({
        "vaf": df['vaf'].values,
        "frag_stat": df['frag_stat'].values,
        "fetal_00": fetal_00,
        "fetal_01": fetal_01,
        "fetal_11": fetal_11,
        "maternal_00": maternal_00,
        "maternal_01": maternal_01,
        "maternal_11": maternal_11
    })

    #homo_alt_idx = df.index[df['vaf'] > 0.975]
    #genotype_df.loc[homo_alt_idx, ['fetal_00', 'fetal_01']] = 0.0
    #genotype_df.loc[homo_alt_idx, 'fetal_11'] = 1.0
    return genotype_df


def main():
    args = parse_args()
    pyro.set_rng_seed(args.seed)

    df = load_and_filter_data(args.input, args.seed)
    f_init = estimate_initial_ff(df)
    frag_means_np = initialize_frag_means(df, f_init)
    frag_means = torch.tensor(frag_means_np, dtype=torch.float32)

    X = torch.tensor(df[['vaf', 'frag_stat']].values, dtype=torch.float32)
    guide = run_inference(X, frag_means, args.nsteps)
    posterior_probs = compute_posterior(X, guide, frag_means)
    posterior_probs=posterior_probs.detach().numpy()
    genotype_df = compute_genotype_probs(df, posterior_probs)
    df_merge=pd.merge(df,genotype_df,how='left',on=["vaf","frag_stat"])
    df_merge.to_csv(args.output, index=False , sep="\t")
    print(f"Genotype probabilities written to: {args.output}")


if __name__ == "__main__":
    main()



