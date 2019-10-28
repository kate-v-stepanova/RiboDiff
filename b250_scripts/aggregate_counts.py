import pandas as pd
import os
import sys
import random


def split(number):
    if number == 0:
        return 0
    if number >= 50:
        min_val = int(number * 0.48)
        max_val = int(number * 0.52)
        a = random.randint(min_val, max_val)
    else:
        a = random.randint(0, number)
    return a   

rna_id = sys.argv[1]
rp_id = sys.argv[2]

sample1 = sys.argv[3]
sample2 = sys.argv[4] # control

# whatever is in the 5th argument, we include non-coding genes too
coding_only = True
noncoding = ['non_coding', 'non-coding', 'noncoding']
if len(sys.argv) >= 6:
    if sys.argv[5] in noncoding :
        coding_only = False
prefix = "coding" if coding_only else "all"

outfile_rna = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ribo_diff/reads/{}_vs_{}_{}.tsv".format(rna_id, sample1, sample2, prefix)
outfile_rp = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ribo_diff/reads/{}_vs_{}_{}.tsv".format(rp_id, sample1, sample2, prefix)

indir_rna = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/alignments/reads_per_gene".format(rna_id)
indir_rp = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/alignments/reads_per_gene".format(rp_id)

gene_names = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/gene_names.txt"
coding_genes_file = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/coding_genes.txt"

sno_rnas_file ="/icgc/dkfzlsdf/analysis/OE0532/static/hg19/snoRNAs/snoRNAs.txt"
sno_rnas = pd.read_csv(sno_rnas_file, sep="|", header=None, names=['gene_id', 'gene_name'])


names_df = pd.read_csv(gene_names, sep="|", header=None, names=["gene_id", "gene_name"], usecols=[1,2])
names_df = names_df.drop_duplicates()

with open(coding_genes_file, 'r') as fff:
    coding_genes = fff.read().splitlines()

#all_samples = {'1': '1_16h_21perO2', '2': '2_16h_1perO2', '3': '3_16h_1perO2', '4': '4_16h_1perO2', '5': '5_48h_1perO2', '6': '6_48h_21perO2'}
#samples = [all_samples[sample1], all_samples[sample2]]
samples = [sample1, sample2]

print("Samples: {}".format(",".join(samples)))

# aggregate
df = None
for sample in samples:
    rna_file = os.path.join(indir_rna, "{}_ReadsPerGene.out.tab".format(sample))
    rp_file = os.path.join(indir_rp, "{}_ReadsPerGene.out.tab".format(sample))
    rna_df = pd.read_csv(rna_file, skiprows=4, names=["gene_id", "RNA_{}".format(sample)], usecols=[0,1], sep="\t")
    rp_df = pd.read_csv(rp_file, skiprows=4, names=["gene_id", "RP_{}".format(sample)], usecols=[0,1], sep="\t")
    df1 = pd.merge(rp_df, rna_df, on="gene_id", how="outer")
    if df is None:
        df = df1
    else:
        df = pd.merge(df, df1, on="gene_id", how="outer")

# filter 
# remove if all zeros
cols = list(df.columns)
cols.remove("gene_id")
df = df.loc[(df[cols]>=20).any(1)]
cols = list(filter(lambda x: x.startswith("RP"), cols))
df = df.loc[(df[cols]>=100).any(1)]

cols = list(df.columns)
cols.remove("gene_id")
for col in cols:
   new_col1 = "R1_" + str(col)
   new_col2 = "R2_" + str(col)
   df[new_col1] = df[col].apply(split)
   df[new_col2] = df[col].sub(df[new_col1], axis=0)
   df = df.drop(columns=[col])

if coding_only:
    df = df.loc[df['gene_id'].isin(coding_genes)]

df = pd.merge(df, names_df, on="gene_id", how="left")
# fill na with gene_id
df['gene_name'] = df['gene_name'].fillna(df['gene_id'])
df = df.loc[~df['gene_name'].str.startswith('MT-')]
df = df.loc[~df['gene_name'].str.startswith('HIST')]
df = df.loc[~df['gene_name'].isin(sno_rnas['gene_name'].tolist())]

df['Entry'] = df['gene_name']
df = df.drop(columns=["gene_name", "gene_id"])  

# reordering
cols = list(df.columns)
cols.remove("Entry")
cols = ["Entry"] + cols
df = df[cols]

print("Writing file: {}".format(outfile_rna))
print("Writing file: {}".format(outfile_rp))
df.to_csv(outfile_rna, sep="\t", header=True, index=False)
df.to_csv(outfile_rp, sep="\t", index=False)


