import os
import sys
import pandas
import glob
import numpy as np

def copy_files():
    files=[]
    for ff in glob.glob('../../HL*'):
        sample = os.path.split(ff)[1]
        cns_file = os.path.join(os.path.abspath(ff),sample + '.readyForPlotting.cns')
        os.system('ln ' + cns_file + ' ' + sample + '.cns')

copy_files()

existing_files = glob.glob('HL*')
existing_files.sort()
samples = map(lambda x: x.split('.')[0],existing_files)

sample_data = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
sample_data = sample_data[sample_data['sample_ID'].isin(samples)]
sample_brca1flfl = sample_data[sample_data['brca1']=='flox/flox']['sample_ID'].tolist()
sample_brca1wt = sample_data[sample_data['brca1']=='WT']['sample_ID'].tolist()

sample_brca1flfl.sort()
sample_brca1wt.sort()

os.system('cnvkit.py heatmap -d ' + ' '.join([i + '.cns' for i in sample_brca1flfl] + [i + '.cns' for i in sample_brca1wt]) + ' -o heatmap.png -x f')

os.system('cnvkit.py heatmap -d ' + ' '.join([i + '.cns' for i in sample_brca1flfl]) + ' -o brca1-flfl-heatmap.png -x f')

os.system('cnvkit.py heatmap -d ' + ' '.join([i + '.cns' for i in sample_brca1wt]) + ' -o brca1-wt-heatmap.png -x f')
