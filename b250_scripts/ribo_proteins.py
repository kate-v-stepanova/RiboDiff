import pandas as pd
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt

project_id = sys.argv[1]

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ribo_diff/te".format(project_id)
ribo_proteins = "/icgc/dkfzlsdf/analysis/OE0532/static/hg19/ribosomal_proteins.txt"
r_prots = pd.read_csv(ribo_proteins, sep="\t", header=None, names=["geneIDs"])
outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/figures/ribo_diff".format(project_id)

df1 = None
for f in glob.glob(os.path.join(indir, "coding*.txt")):
    df = pd.read_csv(f, sep="\t")
    rp = pd.merge(df, r_prots, on="geneIDs", how="inner")
    te_cols = []
    cols = []
    for col in rp.columns:
        if col.startswith('TES'):
            col1 = col.replace('TES', 'TE(').replace('_ctrl', '') + ')'
            cols.append(col1)
            te_cols.append(col)
    te_cols.reverse()
    cols.reverse()
   
    rp = rp[['geneIDs'] + te_cols]
    rp.columns = ['geneIDs'] + cols
    fig = plt.figure()
    plt.title(os.path.basename(f).replace('.txt', '').replace('coding_', ''))
    boxplot = rp.boxplot(column=cols)
    y1 = rp[cols[0]]
    x1 = np.random.normal(1, 0.04, size=len(y1))
    plt.plot(x1, y1, 'r.', alpha=0.2)

    y2 = rp[cols[1]]
    x2 = np.random.normal(2, 0.04, size=len(y2))
    plt.plot(x2, y2, 'r.', alpha=0.2)

    plotname = os.path.join(outdir, 'boxplot.' + os.path.basename(f).replace('.txt', '.pdf'))
    fig.savefig(plotname, format='pdf')
    filename = os.path.join(indir, 'ribosomal_proteins.' + os.path.basename(f))
    rp = rp.round(3)
    rp.to_csv(filename, sep="\t", index=False)

