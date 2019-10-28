#!/usr/bin/env python 
"""
Output result files.
"""

import numpy as np
import cPickle as pickle

def write_result(data, opts):
    """
    Output the result to plain text format.

    @args data: Store all input data and results
    @type data: Class object
    @args opts: Input argument to the main TE function 
    @type opts: Instance
    """
    
    geneIDs = data.geneIDs
    pval = data.pval.astype(str)
    padj = data.padj.astype(str)
    TEctl = data.TEctl.astype(str)
    TEtrt = data.TEtrt.astype(str)
    nameCondA = data.nameCondA
    nameCondB = data.nameCondB
    logFoldChangeTE = data.logFoldChangeTE.astype(str)
    # get counts for RP
    cntRiboNorm = data.countRibo / data.libSizesRibo
    cntRiboMean = np.mean(cntRiboNorm, axis=1)
    log2Ribo = np.log2(cntRiboMean)
    cntRiboMean = cntRiboMean.reshape(len(cntRiboMean), 1)
    log2Ribo = log2Ribo.reshape(len(log2Ribo), 1)
    # get counts for RNA
    cntRnaNorm  = data.countRna  / data.libSizesRna
    cntRnaMean = np.mean(cntRnaNorm, axis=1)
    log2Rna = np.log2(cntRnaMean)
    cntRnaMean = cntRnaMean.reshape(len(cntRnaMean), 1)
    log2Rna = log2Rna.reshape(len(log2Rna), 1)
    meanRnaTrt = data.meanRnaTrt.reshape(len(data.meanRnaTrt), 1)
    meanRnaCtl = data.meanRnaCtl.reshape(len(data.meanRnaCtl), 1)
    meanRiboTrt = data.meanRiboTrt.reshape(len(data.meanRiboTrt), 1)
    meanRiboCtl = data.meanRiboCtl.reshape(len(data.meanRiboCtl), 1)
    if opts.dispDiff:
        dispAdjRibo = data.dispAdjRibo.astype(str)
        dispAdjRna  = data.dispAdjRna.astype(str)
        outNdarrayUnsorted = np.hstack([geneIDs, dispAdjRibo, dispAdjRna, pval, padj, TEctl, TEtrt, logFoldChangeTE, cntRnaMean, cntRiboMean, log2Rna, log2Ribo, meanRnaTrt, meanRnaCtl, meanRiboTrt, meanRiboCtl])
        header = 'geneIDs\tdisperRibo\tdisperRNA\tpval\tpadj\tTE%s\tTE%s\tlog2FC_TE(%s vs %s)\tcntRnaMean\tcntRiboMean\tlog2Rna\tlog2Ribo\tcntNormRna_%s\tcntNormRna_%s\tcntNormRibo_%s\tcntNormRibo_%s' % (nameCondA, nameCondB, nameCondB, nameCondA, nameCondB, nameCondA, nameCondB, nameCondA)
    else:
        dispAdj = data.dispAdj.astype(str)
        outNdarrayUnsorted = np.hstack([geneIDs, dispAdj, pval, padj, TEctl, TEtrt, logFoldChangeTE, cntRnaMean, cntRiboMean, log2Rna, log2Ribo, meanRnaTrt, meanRnaCtl, meanRiboTrt, meanRiboCtl])
        header = 'geneIDs\tdisper\tpval\tpadj\tTE%s\tTE%s\tlog2FC_TE(%s vs %s)\tcntRnaMean\tcntRiboMean\tlog2Rna\tlog2Ribo\tcntNormRna_%s\tcntNormRna_%s\tcntNormRibo_%s\tcntNormRibo_%s' % (nameCondA, nameCondB, nameCondB, nameCondA, nameCondB, nameCondA, nameCondB, nameCondA)

    if opts.rankResult == 0:
        outNdarray = outNdarrayUnsorted.copy()
    else:
        if opts.rankResult == 1:
            idx = np.argsort(data.padj, axis=None)
        elif opts.rankResult == 2:
            idx = np.argsort(data.logFoldChangeTE, axis=None)
        elif opts.rankResult == 3:
            idx = np.argsort(data.geneIDs, axis=None)
        else:
            pass

        outNdarray = outNdarrayUnsorted[idx]
    np.savetxt(opts.outFile, outNdarray, fmt='%s', delimiter='\t', header=header, comments='')

def save_data(data, opts):
    """
    Saving in python object form.

    @args data: Store all input data and results
    @type data: Class object
    @args opts: Input argument to the main TE function 
    @type opts: Instance
    """

    if '.' not in opts.outFile or opts.outFile.endswith('.pkl'):
        pklFile = opts.outFile + '.pkl'
    else:
        pos = opts.outFile.rfind('.')
        if opts.outFile[pos:pos+2] == './':
            outputNamePrefix = opts.outFile
        else:
            outputNamePrefix = opts.outFile[:pos]
        pklFile = outputNamePrefix + '.pkl'
    with open(pklFile, 'wb') as FileOut:
        pickle.dump(data, FileOut, pickle.HIGHEST_PROTOCOL)

