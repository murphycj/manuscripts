import vcf
import sys
import os
import argparse
import pandas
from scipy import stats

def main(args):

    support = pandas.read_csv(args.support,index_col=0)
    coverage = pandas.read_csv(args.coverage,index_col=0)

    samples = pandas.read_excel(args.samples,'samples')
    rnaseq = samples[(samples['sequencing_type']=='RNAseq')]['sample_ID'].tolist()
    rnaseq = list(set(rnaseq).intersection(set(support.columns.tolist())))
    wex = samples[(samples['sequencing_type']=='WEX')]['sample_ID'].tolist()
    wex = list(set(wex).intersection(set(support.columns.tolist())))

    CONTROLS = samples[(samples['tissue']!='primary') & (samples['tissue']!='implant')]['sample_ID'].tolist()
    samples = samples[(samples['tissue']=='primary') | (samples['tissue']=='implant')]
    samples.index = samples['sample_ID']

    primary_tumors = list(set(samples['mouse'].tolist()))
    sample_map = dict(zip(samples['sample_ID'].tolist(),samples['mouse'].tolist()))

    vcf_in = vcf.Reader(open(args.vcf,'r'))
    vcf_out = vcf.Writer(open(args.out,'w'), vcf_in)
    vcf_out_removed = vcf.Writer(open(args.out.replace('filtered','removed'),'w'), vcf_in)

    for v in vcf_in:

        # determine how many of the WEX and RNA-seq samples are
        # mutated and not mutated

        name = str(v.CHROM) + '-' + str(v.POS)

        rnaseq_mut=(support.loc[name,rnaseq]>=3).sum()
        rnaseq_nonmut=(support.loc[name,rnaseq]<3).sum()
        wex_mut=(support.loc[name,wex]>=3).sum()
        wex_nonmut=(support.loc[name,wex]<3).sum()
        pval = stats.fisher_exact([[rnaseq_nonmut,wex_nonmut],[rnaseq_mut,wex_mut]])[1]

        #filter by control

        num_with_min_coverage = 0
        filter_because_control = False

        for s in CONTROLS:
            if s in coverage.columns:
                if args.Cc!=-1 and coverage.ix[name][s] >= args.Cc:
                    num_with_min_coverage+=1

                if args.Cr!=-1 and support.ix[name][s] >= args.Cr:
                    filter_because_control=True

        if args.Cn!=-1 and num_with_min_coverage < args.Cn:
            filter_because_control = True

        #filter by tumor

        filter_because_tumor = False

        primary_tumor_count = {i:False for i in primary_tumors}
        mutant_samples=[]

        for s in v.samples:

            if s.sample not in CONTROLS:

                if s.called:
                    primary_tumor_count[sample_map[s.sample]]=True
                    mutant_samples.append(s.data.AD)

                if support.ix[name][s.sample]>=args.Tr:
                    primary_tumor_count[sample_map[s.sample]]=True

        n_primary = sum(map(lambda x: 1 if x else 0, primary_tumor_count.values()))

        if n_primary >= args.Tn:
            filter_because_tumor=True

        if not filter_because_tumor and not filter_because_control and pval >= 0.01:
            vcf_out.write_record(v)
        else:
            vcf_out_removed.write_record(v)

    vcf_out.close()

parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf',type=str,help='VCF file',required=True)
parser.add_argument('--samples',type=str,help='samples',required=True)
parser.add_argument('--Tn',type=int,help='',required=True)
parser.add_argument('--Tr',type=int,help='',required=True)
parser.add_argument('--Cn',type=int,help='',required=True)
parser.add_argument('--Cc',type=int,help='',required=True)
parser.add_argument('--Cr',type=int,help='',required=True)
parser.add_argument('--coverage',type=str,help='',required=True)
parser.add_argument('--support',type=str,help='',required=True)
parser.add_argument('--out',type=str,help='out vcf file',required=False)
args = parser.parse_args()

main(args=args)
