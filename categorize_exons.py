"""
Usage: python3 categorize_exons.py <reference gff file> <significant exons file (tsv/csv)> <output filename (csv)>

"""

import os
import sys
import pandas as pd
from gtfparse import read_gtf

def categorize_exons(gff, gene, exonID):
    """Categorize significant exons by exon usage (TSS, TTS, Known ES, Novel ES). Each exon falls into 1 category, with the exception that some exons are both TSS and TTS.

    Args:
        gff (df): reference gff file imported as pandas df
        gene (string): Ensembl gene id. 'groupID' in sig exons df (filtered dxr1 from DEXseq results)
        exonID (string): Exon number. 'featureID' in sig exons df (filtered dxr1 from DEXseq results)

    Results:
        final_df (df): row - gene:exonID; columns - TSS, TTS, KES, NES (0 or 1 represents whether an exon falls into each category) 
    """
    exon = exonID.replace('E', '')
    gene_df = gff_data.loc[gff_data['gene_id'].str.contains(gene)]
    # list of unique transcripts per gene
    tscripts = gene_df.loc[gene_df['feature'] == 'exonic_part']['transcripts'].tolist()
    uniq_tscipts = sorted(set([i for sublist in [x.split('+') for x in tscripts] for i in sublist]))
    # rebuild transcripts
    tscript_sets = {}
    for t in uniq_tscipts:
        t_list = gene_df.loc[gene_df['transcripts'].str.contains(t)]['exonic_part_number'].tolist()
        tscript_sets[t] = t_list
    # subset to include entries with exon
    exon_set = {k:v for k,v in tscript_sets.items() if exon in v}
    # categorize exons
    TSS,TTS,KES,NES = (0,0,0,0)
    # 1: TSS [the exon is first in at least one transcript]
    TSS = max([1 if exon == v[0] else 0 for k,v in exon_set.items()])
     # 1: TTS [the exon is last in at least one transcript]
    TTS = max([1 if exon == v[len(v)-1] else 0 for k,v in exon_set.items()])
    if TSS == 0 and TTS == 0:
        # 3: Known exon skipping [this exon absent in at least one transcript]
        if len(exon_set) < len(tscript_sets):
            KES = 1
        # 4: Novel exon skipping [this exon is present in all transcripts]
        elif len(exon_set) == len(tscript_sets):
            NES = 1
    # create summary df
    final_df = pd.DataFrame.from_dict({0: {"Exon":f'{gene}:{exonID}', "TSS":TSS, "TTS":TTS, "KES":KES, "NES":NES}}, orient='index')
    return(final_df)

# input args
gff_file, sig_exons, output_fn = sys.argv[1:4]
# read in gff
gff_data = read_gtf(gff_file)
# read in sig exons df (filtered dxr1)
sig_exons = pd.read_csv(sig_exons, sep = '\t').sort_values(by=['groupID', 'featureID'])
# create results df
sig_exons_results = pd.DataFrame()
for row_index,row in sig_exons.iterrows():
    gene,exonID = row[['groupID', 'featureID']]
    results = categorize_exons(gff_data, gene, exonID)
    sig_exons_results = sig_exons_results.append(results)

# how many exons are both TSS and TTS?
s = sig_exons_results.loc[(sig_exons_results['TSS'] == 1) & (sig_exons_results['TTS'] == 1)]
print(len(s))

# write counts to df to plot in R 
exon_sums = sig_exons_results.sum(axis=0).drop('Exon')
exon_sums.to_csv(output_fn, header=False)

