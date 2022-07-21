import argparse
import sys
import os
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
import glob
import subprocess
import matplotlib.pyplot as plt
import pdb

parser = argparse.ArgumentParser(description = '''Analyse the relationship between the TM-score and pDockQ.''')

parser.add_argument('--pDockQ_df', nargs=1, type= str, default=sys.stdin, help = 'Path to pDockQ scores for complexes.')
parser.add_argument('--TMscore_df', nargs=1, type= str, default=sys.stdin, help = 'Path to TMscores against native complex')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to outdir')


##############FUNCTIONS###############
def pdockq_tmscore(merged, outdir):
    '''Analyse the relationship between pdockq and tmscore
    '''

    av_pdockq = []
    nchains = []
    ninterfaces = []
    models = merged.scored_model.unique()
    tmscores = []
    for model in models:
        sel = merged[merged.scored_model==model]
        av_pdockq.append(sel.pDockQ.mean())
        nchains.append(np.unique(np.concatenate([sel.Chain1.values,sel.Chain2.values])).shape[0])
        ninterfaces.append(len(sel))
        tmscores.append(sel['TM-score'].values[0])

    #Plot relationships
    fig,ax = plt.subplots(figsize=(9/2.54,9/2.54))
    plt.scatter(av_pdockq,tmscores,s=1)
    plt.xlabel('pDockQ')
    plt.ylabel('TM-score')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'pDockQ_vs_TMscore.png',dpi=300)
    plt.close()


    fig,ax = plt.subplots(figsize=(9/2.54,9/2.54))
    plt.scatter(nchains,tmscores,s=1)
    plt.xlabel('# chains')
    plt.ylabel('TM-score')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'nchains_vs_TMscore.png',dpi=300)
    plt.close()

    fig,ax = plt.subplots(figsize=(9/2.54,9/2.54))
    plt.scatter(ninterfaces,tmscores,s=1)
    plt.xlabel('# interfaces')
    plt.ylabel('TM-score')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(outdir+'ninterfaces_vs_TMscore.png',dpi=300)
    plt.close()

    #df
    av_pdockq_df = pd.DataFrame()
    av_pdockq_df['pDockQ']=av_pdockq
    av_pdockq_df['model']=models
    av_pdockq_df = av_pdockq_df.sort_values(by='pDockQ',ascending=False).reset_index()
    top_model = av_pdockq_df.loc[0]
    print('Top pDockQ:',top_model.pDockQ, 'for model',top_model.model)

    pdb.set_trace()

################MAIN###############
#Parse args
args = parser.parse_args()
pDockQ_df = pd.read_csv(args.pDockQ_df[0])
TMscore_df = pd.read_csv(args.TMscore_df[0])
outdir = args.outdir[0]
#Merge dfs
merged = pd.merge(pDockQ_df,TMscore_df,on='scored_model')
pdockq_tmscore(merged, outdir)
#Get best model

pdb.set_trace()
