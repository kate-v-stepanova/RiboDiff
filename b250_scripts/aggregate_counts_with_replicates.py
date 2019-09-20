import pandas as pd
import os
import sys
import random
from functools import reduce

rna_id = sys.argv[1]
rp_id = sys.argv[2]

sample1 = sys.argv[3] # e.g: MDA231_shFTO
sample2 = sys.argv[4] # control, e.g. MDA231_SCR

s1_r1 = sample1 + "1"
s1_r2 = sample1 + "2"
s2_r1 = sample2 + "1"
s2_r2 = sample2 + "2"
samples = [sample1, sample2]

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


print("Samples: {}".format(",".join(samples)))

# aggregate
df = None
for sample in samples:
    rna_file1 = os.path.join(indir_rna, "{}1_ReadsPerGene.out.tab".format(sample))
    rna_file2 = os.path.join(indir_rna, "{}2_ReadsPerGene.out.tab".format(sample))
    rp_file1 = os.path.join(indir_rp, "{}1_ReadsPerGene.out.tab".format(sample))
    rp_file2 = os.path.join(indir_rp, "{}2_ReadsPerGene.out.tab".format(sample))
    rna_df1 = pd.read_csv(rna_file1, skiprows=4, names=["gene_id", "RNA_{}1".format(sample)], usecols=[0,1], sep="\t")
    rna_df2 = pd.read_csv(rna_file2, skiprows=4, names=["gene_id", "RNA_{}2".format(sample)], usecols=[0,1], sep="\t")
    rp_df1 = pd.read_csv(rp_file1, skiprows=4, names=["gene_id", "RP_{}1".format(sample)], usecols=[0,1], sep="\t")
    rp_df2 = pd.read_csv(rp_file2, skiprows=4, names=["gene_id", "RP_{}2".format(sample)], usecols=[0,1], sep="\t")
    dfs = [rna_df1, rp_df1, rna_df2, rp_df2]
    df1 = reduce(lambda left,right: pd.merge(left,right,on='gene_id', how="outer"), dfs)
#    df1 = pd.merge(rp_df, rna_df, on="gene_id", how="outer")
    if df is None:
        df = df1
    else:
        df = pd.merge(df, df1, on="gene_id", how="outer")

df = df.fillna(0)

# filter 
# remove if all zeros
cols = list(df.columns)
cols.remove("gene_id")
df = df.loc[(df[cols]>=20).any(1)]
cols = list(filter(lambda x: x.startswith("RP"), cols))
df = df.loc[(df[cols]>=100).any(1)]

if coding_only:
    df = df.loc[df['gene_id'].isin(coding_genes)]

df = pd.merge(df, names_df, on="gene_id", how="left")
# fill na with gene_id
df['gene_name'] = df['gene_name'].fillna(df['gene_id'])

# remove MT- and HIST genes
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


