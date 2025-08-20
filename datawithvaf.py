import gzip
import pandas as pd

vcf_path=input(str("enter vcf_file :")).strip()

output_path=input(str("enter tsv_vcf_file output:")).strip()
data = []
header = []

with gzip.open(vcf_path, 'rt') as f:
    for line in f:
        if line.startswith('##'):
            continue  # ignorer les lignes de meta-info
        if line.startswith('#CHROM'):
            colnames = line.strip().split('\t')[:5]
            continue
        cols = line.strip().split('\t')
        fmt = cols[8].split(':')   # FORMAT
        sample = cols[9].split(':')  # colonne du premier échantillon
        if len(sample) >= 3:  # vérifier qu'on a au moins 3 valeurs
            vaf = sample[2]
        else:
            vaf = None
        data.append(cols[:5] + [vaf])

# On crée un DataFrame
df = pd.DataFrame(data, columns=colnames + ['VAF'])

# Sauvegarde du résultat
df.to_csv(output_path,sep='\t', index=False)
