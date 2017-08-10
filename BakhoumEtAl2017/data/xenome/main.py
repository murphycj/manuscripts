import os
import glob

fout = open('main.sh','w')
for ss in glob.glob('/scratchLocal/chm2059/*_RNAseq'):
    sample = os.path.split(ss)[1].replace("_RNAseq","")
    for fastq in glob.glob(os.path.join(ss,'*fastq')):
        fout.write('mv ' + fastq + ' /athena/elementolab/scratch/chm2059/from_dat02/chm2059/0816_sam/data/xenome/' + sample + "\n")
fout.close()
