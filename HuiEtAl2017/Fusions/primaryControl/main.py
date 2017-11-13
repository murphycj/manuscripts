import pyensembl
import pandas
import numpy as np

def count_fusions():

    samples = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
    samples = samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0)]
    sample_map = dict(zip(samples['sample_name'].tolist(),samples['sample_ID'].tolist()))

    fusions = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/Fusions/primaryControl/fusions-noNormal.noReadthrough.n2.csv")

    fusions_per_sample = {}

    samples = list(set(fusions['sample'].tolist()))
    for sample in samples:
        fusions_per_sample[sample] = []
    for g in fusions.groupby(['sample','Gene_1_symbol.5end_fusion_partner.','Gene_2_symbol.3end_fusion_partner.']):
        fusions_per_sample[g[0][0]].append(g[0][1] + '-' + g[0][2])

    for sample in fusions_per_sample.keys():
        fusions_per_sample[sample] = set(fusions_per_sample[sample])

    return fusions_per_sample

def count_by_primary():

    fusions = count_fusions()

    samples = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
    samples = samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0)]
    samples['sample_name'] = samples['sample_name'].astype(str)
    samples.index = samples['sample_ID']
    samples = samples[(samples['tissue']!='liver') & (samples['tissue']!='normal breast') & (samples['tissue']!='tail')]
    samples = samples[samples['tissue']=='primary']

    fout = open('fusions-by-primary.tsv','w')
    fout.write('Primary_tumor_number\tSamples\tBRCA1\tTotal_fusions\n')
    for group in samples.groupby(['mouse','primary_tumor_index']):
        if str(group[0][0]) in ['517','1783','1867','2029']:
            primary_tumor_number = str(group[0][0]) + '-' + str(group[0][1])
        else:
            primary_tumor_number = str(group[0][0])
        tmp_samples = ', '.join((group[1]['sample_ID'] + ' (' + group[1]['sequencing_type'] + ')').tolist())

        primary_fusions = set()
        for ss in group[1]['sample_ID']:
            if ss in fusions:
                primary_fusions = primary_fusions.union(fusions[ss])


        fout.write(
            '%s\t%s\t%s\t%s\n' %
            (
                primary_tumor_number,
                tmp_samples,
                group[1]['brca1'].tolist()[0],
                len(primary_fusions)
            )
        )
    fout.close()

count_by_primary()
