import pandas as pd 

data= pd.read_csv("/mnt/B/04_GenMol/04_Logiels/Linux/NGS_nour/dpni_reads/Cardiorun30/extra_Cardiorun30.tsv", sep="\t")
data=data.shift(periods=1)
data.iloc[0,:]=data.columns
data.columns=("Chr","Pos" ,"Ref","Alt","AF","DP","AC")

for col in ["AF","DP","AC"]:
         data[col] = pd.to_numeric(data[col], errors='coerce')
data["AF"] = data.apply(lambda x: x["AC"] / x["DP"] if x["DP"] != 0 else 0, axis=1)
df=pd.read_csv("/mnt/B/04_GenMol/04_Logiels/Linux/NGS_nour/dpni_reads/Cardiorun30/insert_size_ranksum.cleaned.tsv", sep="\t")
#split_df = df["Variant_ID"].str.extract(r'(?P<Chr>[^:]+):(?P<pos>\d+):(?P<Ref>[^>]+)>(?P<Alt>.+)')

# Convert positions to integer
#split_df["pos"] = split_df["pos"].astype(int)
data["Pos"] = data["Pos"].astype(int)
# Concatenate with original dataframe if needed
#df = pd.concat([df, split_df], axis=1)
#df.pop("Variant_ID")
#df.to_csv("/mnt/B/04_GenMol/04_Logiels/Linux/NGS_nour/dpni_reads/variants_split.tsv", sep="\t", index=False)

#print("✅ Output saved to 'variants_split.tsv'")


df_merge=pd.merge(df,data,how='left',on=["Chr","Pos" ,"Ref","Alt"])
DP_data=df_merge.iloc[:,3:]
DP_data.rename(columns={'AC':'alt_count','AF':'vaf','Z_score':'frag_stat','DP':'depth'},inplace=True)
DP_data['frag_stat']=df_merge["Z_score"]
DP_data.to_csv("/mnt/B/04_GenMol/04_Logiels/Linux/NGS_nour/dpni_reads/Cardiorun30/DP_variantsrun30FF.tsv",sep="\t", index=False)



