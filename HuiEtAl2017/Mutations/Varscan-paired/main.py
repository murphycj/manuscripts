import pandas

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
