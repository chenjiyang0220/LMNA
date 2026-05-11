from pomegranate import *
import pandas as pd
import json
import os
import sys
from pybedtools import BedTool
import numpy as np


# paths for LADs
bedgraphs_dir = './bulkCUT/final_20kb_bdg/WT'
outdir = './results/LMNA/WT'


# organization df (org_df) is a file describing the samples
# required columns are cell_type, replicate, chip
# be a tsv/csv etc if you change this command
org_df = pd.read_table('./bulkCUT/LAD_calling/results/sample_info_LA.csv', sep='\t')

# states dict for LADs
states_dict = {
        2:['nonLAD','LAD'],
        3:['nonLAD','T2-LAD','T1-LAD'],
        4:['nonLAD','T3-LAD','T2-LAD','T1-LAD'],
        5:['nonLAD','T4-LAD','T3-LAD','T2-LAD','T1-LAD']
    }
    
encode_blocklist = BedTool('/media/ggj/ggj/CJY/tools/mm10/mm10.blacklist.sorted.bed')


def parse_bedgraph(f, binsize, blocklistbed):
    """
    parses bedgraph file for use in the HMM training
    f: path to bedGraph file
    binsize: resolution of bedgraph file (we used 20kb)
    blocklistbed: BED file of regions to remove, e.g. ENCODE blocklist
    """
    
    # filter non-canonical chromosomes, this will have to be changed
    # if this is applied to nonhuman samples
    canonical_chroms = [ 'chr' + str(x) for x in range(1,23) ] + ['chrY', 'chrX']
    
    # load bedGraph file and convert to BED format
    indf_pre = pd.read_csv(f, sep='\t', low_memory=False, header=None, 
        names=['chrom','start','stop','score']).dropna().query('chrom in @canonical_chroms')
    indf_pre['id'] = indf_pre.index
    
    inbed = BedTool.from_dataframe(indf_pre[['chrom','start','stop','id','score']]).sort()
    inbed_filt = inbed.subtract(blocklistbed).sort()
    
    indf = inbed_filt.to_dataframe()
    indf.columns = ['chrom','start','stop','id','score']
    indf = indf[['chrom','start','stop','score']].copy()

    indf['size'] = indf['stop'] - indf['start']

    # for some reason deepTools sometimes produces bins > specified resolution
    # and that messes up the HMM, so this fixes that
    indf_split = indf.query('size > @binsize').copy() 
    
    if len(indf_split) > 0: 
        to_add = []
        for ix, row in indf_split.iterrows():
            n_intervals = row['size'] / binsize
            start = row['start']
            chrom = row['chrom']
            score = row[f'score']
            for val in range(int(n_intervals)):
                to_add.append(pd.DataFrame({
                    'chrom':[row['chrom']],
                    'start':[start],
                    'stop':[start + binsize],
                    'score':[score]
                }))
                start = start + binsize

        indf = pd.concat([indf.query('size == @binsize')[['chrom','start','stop',f'score']],
                          pd.concat(to_add)]).query('chrom in @canonical_chroms')
    else:
        indf = indf.query('size == @binsize')[['chrom','start','stop',f'score']].query('chrom in @canonical_chroms')
    return(indf)

def assign_categories(dat, cat_names):
    """
    Assign categories to specific HMM states
    based on the mean ChIP-seq coverage
    in bins assigned to that state.
    
    
    cat_names should be provided based on which should be
    assigned to sequecing coverage level in ascending order.
    """
    n_states = len(dat['hmm_pred'].unique())
    n_cats = len(cat_names)
    if n_states != n_cats:
        print(f'{n_states} in data, but {n_cats} provided. '\
              f'Number of categories must match the number '\
              f'of HMM states.')
    else:
        
#         cols = ['score0','score1'] # each cell type has 2 replicates
        cols = ['score']
        
        states = []
        means = []
        
        for state in dat['hmm_pred'].unique():
            states.append(state)
            means.append(np.mean(dat.query(f'hmm_pred == {state}')[cols].mean()))
           
        # sort means and associated states, ascending
        sorted_means = sorted(means)
        means_to_states = dict(zip(means, states))
        
        sorted_states = []
        
        for val in sorted_means:
            sorted_states.append(means_to_states[val])
    
        names_to_states = dict(zip(sorted_states, cat_names))
        dat['category'] = dat['hmm_pred'].replace(names_to_states)
        return dat



for ix, row in org_df.query('chip == "LMNA"').iterrows():
    
    # predict LADs/KDDs for 2-5 states to compare
    for n_states in [2]:
        
        # train one HMM model per cell type
        ct = row['cell_type']
        dat_rep1 = parse_bedgraph(f'{bedgraphs_dir}/{ct}.LMNA.log2.20kb.bedgraph',20000, encode_blocklist)
#         dat_rep2 = parse_bedgraph(f'{bedgraphs_dir}/{ct}_rep2.sorted.bedgraph',20000, encode_blocklist)
        model =  HiddenMarkovModel() 
        dat_model = model.from_samples(NormalDistribution, n_components=n_states, X=[dat_rep1['score'].tolist()], algorithm='baum-welch', verbose=True, n_jobs=20)
        dat = dat_rep1
        dat['hmm_pred'] = dat_model.predict(dat[['score']].to_numpy(),
                                        algorithm='viterbi')[1:]
        trained_model_json = dat_model.to_json()
        if not os.path.exists(f'{outdir}/{ct}/hmms_{n_states}states'):
            os.makedirs(f'{outdir}/{ct}/hmms_{n_states}states')
            os.makedirs(f'{outdir}/{ct}/hmm_calls_{n_states}states')
            os.makedirs(f'{outdir}/{ct}/BED_files_{n_states}states')
            
        # save the trained model
        with open(f'{outdir}/{ct}/hmms_{n_states}states/{ct}.json', 'w') as outfile:
            json.dump(trained_model_json, outfile)
            
        # save calls per bin
        dat.to_csv(f'{outdir}/{ct}/hmm_calls_{n_states}states/{ct}.tsv',
              sep='\t', index=False)
        
        # assign LAD/KDD categories per bin, and save BED file
        dat = pd.read_table(f'{outdir}/{ct}/hmm_calls_{n_states}states/{ct}.tsv', names=['chrom','start','stop','score','hmm_pred'],
                           header=0) 
        dat_w_cats = assign_categories(dat, states_dict[n_states])
        for cat in dat_w_cats['category'].unique():
            bed = BedTool.from_dataframe(dat_w_cats.query('category == @cat')).sort().merge(d=1)
            bed.to_dataframe().to_csv(f'{outdir}/{ct}/BED_files_{n_states}states/{ct}_{cat}s.bed',
                             sep='\t', header=False, index=False)