params.samtools = '/home/chm2059/chm2059/lib/samtools-1.3.1/samtools'
params.access = '/home/chm2059/chm2059/data/refdata/GRCm38/cnvkit/access-10kb.GRCm38.bed'
params.reference = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'

params.my_targets_sureselect = '/home/chm2059/chm2059/data/refdata/GRCm38/cnvkit/SureSelect/my_targets.bed'
params.my_antitargets_sureselect = '/home/chm2059/chm2059/data/refdata/GRCm38/cnvkit/SureSelect/my_antitargets.bed'

params.my_targets_nimblegen = '/home/chm2059/chm2059/data/refdata/GRCm38/cnvkit/NimbleGen/my_targets.bed'
params.my_antitargets_nimblegen = '/home/chm2059/chm2059/data/refdata/GRCm38/cnvkit/NimbleGen/my_antitargets.bed'

samples = Channel
  .fromFilePairs("/home/chm2059/chm2059/elementoCantley_2014_9_29/results/WEX/BWAGATK/HL*/HL*.gatk.{bam,bam.bai}")
  .filter{
    it[1].size()==2
  }
  .flatMap{
    it -> [[it[0], it[1][0], it[1][1]]]
  }

process cnvkit_coverage {

  storeDir "${baseDir}/${sample}"
  scratch true
  executor 'sge'
  stageInMode 'copy'
  clusterOptions '-l h_vmem=3G -pe smp 4  -l h_rt=72:00:00 -l os=rhel6.3'

  input:
    set sample, bam, bam_bai from samples

  output:
    set sample, file('*.targetcoverage.cnn'), file('*.antitargetcoverage.cnn') into cnvkit_coverage_out

  """

  ${params.samtools} view -h -q 30 \
    ${bam} \
    | awk \' \$7==\"=\" || \$1 ~ /^@/ { print \$0 }\' \
    | ${params.samtools} view -h -bS - > q30diff.bam
  ${params.samtools} index q30diff.bam

  cnvkit.py coverage \
    -p 4 \
    q30diff.bam ${params.my_targets_sureselect} \
    -o ${sample}.targetcoverage.cnn

  cnvkit.py coverage \
    -p 4 \
    q30diff.bam ${params.my_antitargets_sureselect} \
    -o ${sample}.antitargetcoverage.cnn
  """
}
