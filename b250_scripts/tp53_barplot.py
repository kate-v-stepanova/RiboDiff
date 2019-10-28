import pandas as pd
import numpy as np
import glob
import os
import sys
import matplotlib.pyplot as plt

project_id = sys.argv[1]

indir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/ribo_diff/te".format(project_id)
outdir = "/icgc/dkfzlsdf/analysis/OE0532/{}/analysis/output/figures/ribo_diff".format(project_id)

for f in glob.glob(os.path.join(indir, "coding*.txt")):
    df = pd.read_csv(f, sep="\t")
    te_cols = []
    cols = []
    df = df.loc[df['geneIDs'] == 'TP53']
    for col in df.columns:
        if col.startswith('TES'):
            col1 = col.replace('TES', 'TE(').replace('_ctrl', '') + ')'
            cols.append(col1)
            te_cols.append(col)
    te_cols.reverse()
    cols.reverse()
    
    tp53 = pd.DataFrame(columns=['sample','te'])

    tp53['sample'] = cols
    tp53['te'] = [df.iloc[0][te_cols[0]], df.iloc[0][te_cols[1]]]
    tp53.index = tp53['sample'].tolist()
   
    barplot = tp53.plot(kind='bar', width=0.3, rot=0, legend=False)
    plt.title(os.path.basename(f).replace('.txt', '').replace('coding_', ''))


    plotname = os.path.join(outdir, 'tp53.' + os.path.basename(f).replace('.txt', '.pdf').replace('coding_', ''))
    plt.savefig(plotname, format='pdf')
    filename = os.path.join(indir, 'ribosomal_proteins.' + os.path.basename(f))

