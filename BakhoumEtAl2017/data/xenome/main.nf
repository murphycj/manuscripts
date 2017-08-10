#!/usr/bin/env nextflow

bams = Channel.create()

bam1 = Channel
  .fromFilePairs("/home/chm2059/chm2059/0816_sam/data/Samuel4427_2017_03_03/*_L00*_R1_001.fastq.gz", size: 1)
  .map{[it[0].split('_')[0],it[1][0]]}

bam2 = Channel
  .fromFilePairs("/home/chm2059/chm2059/0816_sam/data/Samuel4427_2017_02_08/*_L00*_R1_001.fastq.gz", size: 1)
  .map{[it[0].split('_')[0],it[1][0]]}

bam1.concat(bam2).filter{!(it[0] in['11','12','13','14','15','16','17','18'])}.set{bams}

process xenome {

  input:
    set prefix, file(read1), from bams

  """
  mkdir -p /athena/elementolab/scratch/chm2059/from_dat02/chm2059/0816_sam/data/xenome/${prefix}
  mkdir -p /scratchLocal/chm2059/${prefix}_RNAseq
  cat ${read1} > /scratchLocal/chm2059/${prefix}_RNAseq/R1.fastq.gz

  singularity exec \
    /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/singularity_images/murphycj-Singularity_gossamer-master.img \
    /usr/local/bin/xenome classify \
    -P /scratchLocal/chm2059/GRCm38.GRCh38/GRCm38.GRCh38 \
    -T 1 \
    -i /scratchLocal/chm2059/${prefix}_RNAseq/R1.fastq.gz \
    --output-filename-prefix /scratchLocal/chm2059/${prefix}_RNAseq/${prefix}  > /scratchLocal/chm2059/${prefix}_RNAseq/xenome.log
  """
}
