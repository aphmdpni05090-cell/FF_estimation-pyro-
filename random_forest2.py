import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, roc_curve, precision_score, recall_score, f1_score




# --- Load and clean data ---
data= pd.read_csv("extracted_input.tsv", sep="\t")
data=data.shift(periods=1)
data.iloc[0,:]=data.columns
data.columns=("Chr","Start" ,"Ref","Alt", "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "DP")
df0=pd.read_csv("proband.annovar.hg38_multianno.txt", sep="\t")
df0["AF_popmax"] = pd.to_numeric(df0["AF_popmax"], errors='coerce')
df1=df0[df0["AF_popmax"]>0]
df1=df1[["Chr","Start","Ref","Alt","AF_popmax"]]
df1["labels"]=0
df1["labels"]=(df1["AF_popmax"] > 0.1).astype(int)
df1["Start" ] = pd.to_numeric(df1["Start" ], errors='coerce')
data["Start" ] = pd.to_numeric(data["Start" ], errors='coerce') 
df2=pd.merge(data,df1,how='right',on=["Chr","Start","Ref","Alt"])


# --- Convert and clean data ---
cols_to_numeric = [
    "AF_popmax", "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "BaseQRankSum", "DP"
]
for col in cols_to_numeric:
    df2[col] = pd.to_numeric(df2[col], errors='coerce')


# --- Add Indel binary feature ---
df2["Indel"] = df2.apply(lambda row: 1 if len(str(row["Ref"])) != len(str(row["Alt"])) else 0, axis=1)

# --- Feature selection ---
features = [ "SOR", "MQ", "MQRankSum", "ReadPosRankSum", "DP","AF_popmax", "Indel",]

X = df2[features].values
y = df2["labels"].astype(int).values


# --- Split Data (Stratified) ---
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, stratify=y, random_state=42)

# --- Train Classifier ---
clf = RandomForestClassifier(n_estimators=800, random_state=42, n_jobs=-1)
clf.fit(X_train, y_train)

# --- Evaluation on Test Set ---
y_probs = clf.predict_proba(X_test)[:, 1]
y_preds = clf.predict(X_test)

# Metrics
roc_auc = roc_auc_score(y_test, y_probs)
precision = precision_score(y_test, y_preds)
recall = recall_score(y_test, y_preds)
f1 = f1_score(y_test, y_preds)

print(f"✅ Evaluation Metrics:")
print(f" - ROC AUC: {roc_auc:.4f}")
print(f" - Precision: {precision:.4f}")
print(f" - Recall: {recall:.4f}")
print(f" - F1 Score: {f1:.4f}")



# --- Predict on Full Dataset ---
df2["rf_score"] = clf.predict_proba(X)[:, 1]

# --- Save Scored Data ---
df2.to_csv("scored_variants_with_rf.tsv", sep="\t", index=False)
fpr, tpr, thresholds = roc_curve(y_test, y_probs)

target_sensitivity = 0.99
idx = np.argmax(tpr >= target_sensitivity)
cutoff_score = thresholds[idx]
print(f"📌 Cutoff for {target_sensitivity*100:.0f}% sensitivity: {cutoff_score:.4f}")


df2["FILTER"] = df2["rf_score"].apply(lambda x: "PASS" if x >= cutoff_score else "FAIL")
df2[df2["FILTER"]=="FAIL"].to_csv("failed_scored_variants.tsv", sep="\t", index=False)
df2[df2["FILTER"]=="PASS"].to_csv("filtered_scored_variants.tsv", sep="\t", index=False)
print("📈 Scored variants written to 'scored_variants_with_rf.tsv'")
