params.samtools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/samtools-1.1/samtools'
params.bed_nimblegen = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/Mouse_Exome_Design/110624_MM10_exome_L2R_D02_EZ_HX1-v4.bed'
params.bed_sureselect = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/Mouse_Exome_Design/S0276129_Regions-mm10-v2.bed'

params.copywriter = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CopywriteR/main.R'

initial_brca1_flfl = [
  'HL-WES-01','HL-WES-02','HL-WES-03','HL-WES-04','HL-WES-05',
  'HL-WES-06','HL-WES-07','HL-WES-08','HL-WES-09','HL-WES-10',
  'HL-WES-11','HL-WES-12','HL-WES-13','HL-WES-14','HL-WES-15','HL-WES-16']


bams = Channel
  .fromFilePairs("/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/HL*/HL*.gatk.{bam,bam.bai}")
  .filter{
    it[1].size()==2
  }
  .flatMap{
    it -> [[it[0], it[1][0], it[1][1]]]
  }
  .filter{
    ! (it[0] in initial_brca1_flfl)
  }
  .spread([1,2,4])


process copywriter {

  storeDir "${baseDir}/SD${SD}/${prefix}"
  maxForks 30

  input:
    set prefix, bam, bai, SD from bams

  output:
    set file("${prefix}_result*") into copywriter_out

  """

  Rscript ${params.copywriter} \
    -bam ${bam} \
    -bed ${params.bed_sureselect} \
    -sample ${prefix} \
    -window 40000 \
    -SD ${SD}
  """
}
