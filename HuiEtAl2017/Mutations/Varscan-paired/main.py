import pandas

def count_mutations():

    mutations = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/Mutations/Varscan-paired/somatic.nodbsnp.varscan.filtered.allEffects.primaryControl.csv")

    mutations_per_sample = {}
    nonsilent_mutations_per_sample = {}
    silent_mutations_per_sample = {}

    samples = mutations.columns.tolist()[10:]
    for sample in samples:
        mutations_per_sample[sample] = []
        nonsilent_mutations_per_sample[sample] = []
        silent_mutations_per_sample[sample] = []

    for g in mutations.groupby(['CHROM','POS','REF','ALT']):
        tmp = '-'.join([str(i) for i in g[0]])
        effect = g[1]['EFFECT'].tolist()[0]

        issilent = False

        if effect in ['3_prime_UTR_variant','5_prime_UTR_variant','intergenic_region','intragenic_variant','intron_variant','non_coding_transcript_exon_variant','synonymous_variant','splice_region_variant&synonymous_variant','splice_region_variant&intron_variant','splice_region_variant']:
            issilent=True


        tmp_index = g[1].index[0]
        for ss in samples:
            if mutations.loc[tmp_index,ss]!='-':
                if issilent:
                    silent_mutations_per_sample[ss].append(tmp)
                else:
                    nonsilent_mutations_per_sample[ss].append(tmp)

                mutations_per_sample[ss].append(tmp)
    for sample in mutations_per_sample.keys():
        mutations_per_sample[sample] = set(mutations_per_sample[sample])
        silent_mutations_per_sample[sample] = set(silent_mutations_per_sample[sample])
        nonsilent_mutations_per_sample[sample] = set(nonsilent_mutations_per_sample[sample])

    return mutations_per_sample, nonsilent_mutations_per_sample, silent_mutations_per_sample

def count_by_primary():

    mutations, nonsilent_mutations, silent_mutations = count_mutations()

    samples = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
    samples = samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0)]
    samples['sample_name'] = samples['sample_name'].astype(str)
    samples.index = samples['sample_ID']
    samples = samples[(samples['tissue']!='liver') & (samples['tissue']!='normal breast') & (samples['tissue']!='tail')]
    samples = samples[samples['tissue']=='primary']

    fout = open('mutations-by-primary.tsv','w')
    fout.write('primary_tumor_number\tsamples\tbrca1\ttotal_mutations\tsilent_mutations\tnonsilent_mutations\n')
    for group in samples.groupby(['mouse','primary_tumor_index']):
        if str(group[0][0]) in ['517','1783','1867','2029']:
            primary_tumor_number = str(group[0][0]) + '-' + str(group[0][1])
        else:
            primary_tumor_number = str(group[0][0])
        tmp_samples = ', '.join((group[1]['sample_ID'] + ' (' + group[1]['sequencing_type'] + ')').tolist())

        primary_mutations = set()
        primary_mutations_silent = set()
        primary_mutations_nonsilent = set()

        for ss in group[1]['sample_ID']:
            if ss in mutations:
                primary_mutations = primary_mutations.union(mutations[ss])
                primary_mutations_nonsilent = primary_mutations_nonsilent.union(nonsilent_mutations[ss])
                primary_mutations_silent = primary_mutations_silent.union(silent_mutations[ss])


        fout.write(
            '%s\t%s\t%s\t%s\t%s\t%s\n' %
            (
                primary_tumor_number,
                tmp_samples,
                group[1]['brca1'].tolist()[0],
                len(primary_mutations),
                len(primary_mutations_silent),
                len(primary_mutations_nonsilent)
            )
        )
    fout.close()

def print_pairs():

    samples = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
    samples = samples[samples['sample_ID']!='HL31']
    samples = samples[samples['sample_ID']!='HL159']
    samples = samples[samples['sample_ID']!='HL22']

    tumors = samples[(samples['tissue']=='primary') | (samples['tissue']=='implant')]
    tails = samples[samples['tissue']=='tail']


    tails['mouse'] = tails['mouse'].astype(str)

    fout = open('paired.csv','w')

    for i in tumors.index:
        mouse = str(tumors.loc[i]['mouse'])
        tumor_sample = tumors.loc[i]['sample_ID']
        tumor = ''
        tail = ''
        if mouse in tails['mouse'].tolist():
            controls = tails[tails['mouse']==mouse]
            control_sample = controls.iloc[0]['sample_ID']
            if tumors.loc[i]['sequencing_type'] == 'WEX':
                tumor = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/%s/%s.gatk.bam' % (tumor_sample,tumor_sample)
            else:
                tumor = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/RNAseq/STARGATK/%s/%s.gatk.bam' % (tumor_sample,tumor_sample)

            if control_sample == 'HL-WES-31':
                control = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/%s/%s.gatk.bam' % ('Liver','Liver')
                control_sample = 'Liver'
            else:
                control = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/%s/%s.gatk.bam' % (control_sample,control_sample)
        else:
            if tumors.loc[i]['sequencing_type'] == 'WEX':
                tumor = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/%s/%s.gatk.bam' % (tumor_sample,tumor_sample)
            else:
                tumor = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/RNAseq/STARGATK/%s/%s.gatk.bam' % (tumor_sample,tumor_sample)

            control = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/%s/%s.gatk.bam' % ('Liver','Liver')
            control_sample = 'Liver'

        fout.write(
            '%s,%s,%s,%s.bai,%s,%s.bai\n' %
            (tumor_sample,control_sample,tumor,tumor,control,control)
        )
    fout.close()

count_by_primary()
