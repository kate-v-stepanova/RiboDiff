#!/bin/bash

project_id=$1
rna_id=$project_id
rp_id=$2
te_script="/icgc/dkfzlsdf/analysis/OE0532/software/RiboDiff/scripts/TE.py"
indir="/icgc/dkfzlsdf/analysis/OE0532/$project_id/analysis/output/ribo_diff/reads"
metadir="/icgc/dkfzlsdf/analysis/OE0532/$project_id/analysis/input/ribo_diff"
outdir="/icgc/dkfzlsdf/analysis/OE0532/$project_id/analysis/output/ribo_diff/te"
rp_outdir="/icgc/dkfzlsdf/analysis/OE0532/$rp_id/analysis/output/ribo_diff/te"
plotdir="/icgc/dkfzlsdf/analysis/OE0532/$project_id/analysis/output/figures/ribo_diff"
plotdir_rp="/icgc/dkfzlsdf/analysis/OE0532/$rp_id/analysis/output/figures/ribo_diff"

mkdir -p $outdir
mkdir -p $plotdir
mkdir -p $rp_outdir
mkdir -p $plotdir_rp

for f in $(ls $indir/*_coding.tsv); do 
    fn=$(basename $f); 
    fn=${fn%_coding.tsv}; 
    echo "bsub python $te_script -e $metadir/${fn}.metafile.csv -c $f -o $outdir/coding_${fn}.txt -p 1";
done

for f in $(ls $indir/*_all.tsv); do 
    fn=$(basename $f); 
    fn=${fn%_all.tsv}; 
    echo "bsub python $te_script -e $metadir/${fn}.metafile.csv -c $f -o $outdir/all_${fn}.txt -p 1";
done

echo ""
echo "WAIT untill the previous jobs are done, then run the following:"
echo ""

echo "rm $outdir/*.pkl"
echo "cp $outdir/* $rp_outdir" 
echo "mv $outdir/*.pdf $plotdir"
echo "mv $rp_outdir/*.pdf $plotdir_rp"
echo "python /icgc/dkfzlsdf/analysis/OE0532/software/diricore/utils/fix_excel.py $project_id"
echo "python /icgc/dkfzlsdf/analysis/OE0532/software/diricore/utils/fix_excel.py $rp_id"

