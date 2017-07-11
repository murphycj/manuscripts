params.reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.reference_dict = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.dict'
params.reference_fai = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'
params.samtools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/samtools-1.3.1/samtools'
params.seqpy = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/'
params.bed_file = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/data/Mouse_Exome_Design/110624_S0276129.bed'


rnaseq_bams = Channel.empty()
wex_bams = Channel.empty()
bams = Channel.empty()
test = Channel.empty()

Channel
  .fromFilePairs("/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/RNAseq/STARGATK/HL*/HL*.gatk.{bam,bam.bai}", size: 2)
  .flatMap{
    it -> [[it[0], it[1][0], it[1][1]]]
  }
  .set{
    rnaseq_bams
  }

Channel
  .fromFilePairs("/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/HL*/HL*.gatk.{bam,bam.bai}", size: 2)
  .flatMap{
    it -> [[it[0], it[1][0], it[1][1]]]
  }
  .set{
    wex_bams
  }

test.concat( rnaseq_bams, wex_bams ).set{bams}


process count_context {

  storeDir "${baseDir}/${s}"
  maxForks 45


  input:
    set s, bam, bam_bai from bams

  output:
    file("*.csv") into count_context_out
    val s into count_context_out_sample

  """
  ${params.samtools} mpileup \
    -Q 15 -q 15 -d10000 -f ${params.reference} -l ${params.bed_file} \
    ${bam} \
    | python /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/count_sequenced_kmer.py \
    --n 12 \
    --out ${s}.count.csv
  """
}

process combine_context {

  storeDir "${baseDir}"

  input:
    file counts from count_context_out.toList()
    val s from count_context_out_sample.toList()

  output:
    file('counts.csv') into combine_context_out

  """
  python ${params.seqpy}/aggregate_trinucleotide_counts.py \
    --counts ${counts} \
    --samples ${s} \
    --out counts.csv

  """
}
