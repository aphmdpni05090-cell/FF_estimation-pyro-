import pandas as pd 

vcf_vaf_path=input(str("enter the vcf with vaf file")).strip()
insert_size_test_path=input(str("enter the insert size test path  file")).strip()
output_path=input(str("enter the output path  file name")).strip()
data= pd.read_csv(vcf_vaf_path, sep="\t")
#data=data.shift(periods=1)
#data.iloc[0,:]=data.columns
#data.columns=("Chr","Pos" ,"Ref","Alt","VAF")

data.rename(columns={'#CHROM':'Chr','POS':'Pos','REF':'Ref','ALT':'Alt'},inplace=True)
df=pd.read_csv(insert_size_test_path, sep="\t")
#split_df = df["Variant_ID"].str.extract(r'(?P<Chr>[^:]+):(?P<pos>\d+):(?P<Ref>[^>]+)>(?P<Alt>.+)')
data.drop('ID', axis=1)
# Convert positions to integer
#split_df["pos"] = split_df["pos"].astype(int)
data["Pos"] = data["Pos"].astype(int)
# Concatenate with original dataframe if needed
#df = pd.concat([df, split_df], axis=1)
#df.pop("Variant_ID")
#df.to_csv("/mnt/B/04_GenMol/04_Logiels/Linux/NGS_nour/dpni_reads/variants_split.tsv", sep="\t", index=False)

#print("✅ Output saved to 'variants_split.tsv'")


df_merge=pd.merge(df,data,how='left',on=["Chr","Pos" ,"Ref","Alt"])
DP_data=df_merge.iloc[:,:4]
DP_data["VAF"]=df_merge["VAF"]
DP_data.rename(columns={'VAF':'vaf','Z_score':'frag_stat'},inplace=True)
DP_data['frag_stat']=df_merge["Z_score"]
DP_data.to_csv(output_path,sep="\t", index=False)



