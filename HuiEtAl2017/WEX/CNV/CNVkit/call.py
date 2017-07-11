import os
import sys

import argparse
import pandas
import numpy as np


def main(args):
    #purity = pandas.read_csv(args.purity,index_col=0)
    #purity = purity.ix[args.sample]['purity']

    cns = pandas.read_table(args.cns,sep='\t')
    cnr = pandas.read_table(args.cnr,sep='\t')

    log2s = cnr["log2"].tolist()
    log2s.sort()

    start = int(len(log2s)*0.25)
    end = int(len(log2s)*0.75)

    log2s_sub = log2s[start:end]

    copy_gain_threshold = np.median(log2s_sub) + 2 * np.std(log2s_sub)
    amplification_threshold = np.median(log2s_sub) + 6 * np.std(log2s_sub)
    hetdel_threshold = np.median(log2s_sub) - 2.5 * np.std(log2s_sub)
    homdel_threshold = np.median(log2s_sub) - 7 * np.std(log2s_sub)

    cns_new = os.path.split(args.cns)[1]
    cns_new = cns_new.replace('.cns','.call.cns')

    cnr_new = os.path.split(args.cnr)[1]
    cnr_new = cnr_new.replace('.cnr','.call.cnr')

    os.system(
        'python ~/.local/bin/cnvkit.py call -m threshold -t=%s,%s,%s,%s -g f -o %s %s' %
        (
            homdel_threshold,
            hetdel_threshold,
            copy_gain_threshold,
            amplification_threshold,
            cns_new,
            args.cns
        )
    )
    os.system(
        'python ~/.local/bin/cnvkit.py call -m threshold -t=%s,%s,%s,%s -g f -o %s %s' %
        (
            homdel_threshold,
            hetdel_threshold,
            copy_gain_threshold,
            amplification_threshold,
            cnr_new,
            args.cnr
        )
    )

    cns = pandas.read_table(cns_new,sep='\t')
    cns.loc[cns['log2']>0,'log2'] = cns.loc[cns['log2']>0,'log2']/amplification_threshold
    cns.loc[cns['log2']<0,'log2'] = cns.loc[cns['log2']<0,'log2']/(-1*homdel_threshold)
    os.system('rm ' + cns_new)
    cns.to_csv(cns_new,sep='\t',index=None)

    cns["chromosome"] = cns["chromosome"].apply(str)
    cns = cns[cns["chromosome"]!="Y"]
    cns.loc[cns['cn']==2,'log2']=0
    tmp = cns_new.replace('.call','').replace('.corrected','').replace('.cns','.readyForPlotting.cns')
    cns.to_csv(tmp,sep='\t',index=None)

    cnr = pandas.read_table(cnr_new,sep='\t')
    cnr.loc[cnr['log2']>0,'log2'] = cnr.loc[cnr['log2']>0,'log2']/amplification_threshold
    cnr.loc[cnr['log2']<0,'log2'] = cnr.loc[cnr['log2']<0,'log2']/(-1*homdel_threshold)
    os.system('rm ' + cnr_new)
    cnr.to_csv(cnr_new,sep='\t',index=None)

    cnr["chromosome"] = cnr["chromosome"].apply(str)
    cnr = cnr[cnr["chromosome"]!="Y"]
    tmp = cnr_new.replace('.call','').replace('.corrected','').replace('.cnr','.readyForPlotting.cnr')
    cnr.to_csv(tmp,sep='\t',index=None)


parser = argparse.ArgumentParser(description='')
parser.add_argument('--sample',type=str,help='Sample name',required=True)
parser.add_argument('--cns',type=str,help='.cns file',required=True)
parser.add_argument('--cnr',type=str,help='.cnr file',required=True)
#parser.add_argument('--purity',type=str,help='File containing purity estimates',required=True)
args = parser.parse_args()

main(args=args)
