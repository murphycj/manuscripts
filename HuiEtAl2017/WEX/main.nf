#!/usr/bin/env nextflow

initial_brca1_flfl = [
  'HL-WES-01','HL-WES-02','HL-WES-03','HL-WES-04','HL-WES-05',
  'HL-WES-06','HL-WES-07','HL-WES-08','HL-WES-09','HL-WES-10',
  'HL-WES-11','HL-WES-12','HL-WES-13','HL-WES-14','HL-WES-15','HL-WES-16']

params.fastqs = "/home/chm2059/chm2059/elementoCantley_2014_9_29/data/WEX/HL*"

params.capture_file = '/home/chm2059/chm2059/data/Mouse_Exome_Design/110624_MM10_exome_L2R_D02_EZ_HX1-v3.bed'
params.capture_file_sureselect = '/home/chm2059/chm2059/data/Mouse_Exome_Design/S0276129_Regions-mm10-v2.bed'

params.bedtools_genome = '/home/chm2059/chm2059/data/refdata/GRCm38/bedtools_windows/GRCm38.txt'
params.bedtools_intervals = '/home/chm2059/chm2059/data/refdata/GRCm38/bedtools_windows/windows.10000.bed'
params.bedtools_intervals_trp53 = '/home/chm2059/chm2059/elementoCantley_2014_9_29/data/Trp53.floxed.bed'

params.reference = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.reference_dict = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.dict'
params.reference_fai = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'

params.bwa = '/home/chm2059/chm2059/lib/bwa-0.7.12/bwa'
params.dbsnp = '/home/chm2059/chm2059/data/refdata/mm10/mgp.v4.snps.dbSNP.vcf'
params.dexseq_gff = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.dexseq.gff'
params.gtf = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.gtf'

params.pindel = '~/chm2059/lib/pindel/pindel'
params.java = '/home/chm2059/chm2059/lib/jre1.8.0_25/bin/java'
params.varscan = '/home/chm2059/chm2059/lib/VarScan.v2.3.9.jar'
params.fastqc = '/home/chm2059/chm2059/lib/FastQC/fastqc'
params.aggregate_fastqc = '/home/chm2059/chm2059/lib/seqpy/bin/aggregate_fastqc.py'
params.seqpy = '/home/chm2059/chm2059/lib/seqpy/bin/'
params.samtools = '/home/chm2059/chm2059/lib/samtools-1.1/samtools'
params.picard = '~/chm2059/lib/picard-tools-1.137/picard.jar'
params.gatk = '~/chm2059/lib/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
params.picard = '/home/chm2059/chm2059/lib/picard-tools-1.137/picard.jar'
params.flagstat = '~/chm2059/lib/seqpy/bin/flagstat.py'

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()


Channel
  .fromPath(params.fastqs, type: 'dir')
  .map { path ->
       (prefix, R1_files, R2_files) = getFastqFiles(path)
       tuple(prefix, R1_files, R2_files)
  }
  .separate(fqfiles, fqfiles2, fqfiles3){it->[it,it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 10

  input:
    set prefix, file(read1), file(read2) from fqfiles

  output:
    file('*zip') into fastqc_results
    val prefix into samples

  """
  zcat ${read1} ${read2} | ${params.fastqc} stdin --outdir=.
  mv stdin_fastqc.zip ${prefix}.zip
  """
}

process combine_fastqc_results {

  storeDir "${baseDir}/fastqc/"

  input:
    file results from fastqc_results.toList()
    val all_samples from samples.toList()

  output:
    file('HL_WEX_Fastqc.csv') into fastqc_combine_output

  """
  python ${params.aggregate_fastqc} --out HL_WEX_Fastqc.csv --files ${results} --samples ${all_samples.join(" ")}
  """
}

process process_gatk {
  scratch true
  executor 'sge'
  stageInMode 'copy'
  clusterOptions '-l h_vmem=8G -pe smp 4 -l h_rt=96:00:00 -l os=rhel5.4|rhel6.3'

  storeDir "${baseDir}/BWAGATK/${prefix}"

  input:
    set prefix, file(read1), file(read2) from fqfiles3

  output:
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out2
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_insert_size
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_flagstat
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_gatk_coverage
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_varscan_snp
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_varscan_indel
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_dexseq_count
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_bedtools_count
    file('*.gatk.bam') into gatk_out_gatk_coverage_bam
    file('*.gatk.bam.bai') into gatk_out_gatk_coverage_bai

  """
  cat ${read1} > ${prefix}.R1.fastq.gz
  cat ${read2} > ${prefix}.R2.fastq.gz
  gunzip ${prefix}.R1.fastq.gz
  gunzip ${prefix}.R2.fastq.gz
  ${params.bwa} aln -t 4 ${params.reference} ${prefix}.R1.fastq > R1.sai
  ${params.bwa} aln -t 4 ${params.reference} ${prefix}.R2.fastq > R2.sai
  ${params.bwa} sampe ${params.reference} R1.sai R2.sai ${prefix}.R1.fastq ${prefix}.R2.fastq > aln-pe.sam
  ${params.samtools} view -Sb aln-pe.sam > tmp.bam
  ${params.samtools} sort tmp.bam sorted
  ${params.samtools} index sorted.bam

  rm tmp.bam
  rm aln-pe.sam

  ${params.java} -Xmx10g -jar \
    ${params.picard} AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=LENIENT \
    I=sorted.bam \
    O=tmp.bam \
    ID=${prefix} \
    LB=1 \
    PL=illumina \
    PU=1 \
    RGSM=test

  ${params.samtools} index tmp.bam

  ${params.java} -Xmx10g -jar \
    ${params.picard} MarkDuplicates \
    VALIDATION_STRINGENCY=LENIENT \
    I=tmp.bam \
    O=tmp2.bam \
    REMOVE_DUPLICATES=true \
    M=duplicates.bam

  ${params.samtools} index tmp2.bam
  rm tmp.bam*

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T RealignerTargetCreator \
    -R ${params.reference} \
    -I tmp2.bam \
    -U ALLOW_N_CIGAR_READS \
    -o forIndelRealigner.intervals \
    -nt 4

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T IndelRealigner \
    -R ${params.reference} \
    -I tmp2.bam \
    -targetIntervals forIndelRealigner.intervals \
    -o tmp.bam

  ${params.samtools} index tmp.bam

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T BaseRecalibrator \
    -R ${params.reference} \
    -I tmp.bam \
    -o recal_data.table \
    -knownSites ${params.dbsnp} \
    -nct 4

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T PrintReads \
    -R ${params.reference} \
    -I tmp.bam \
    --BQSR recal_data.table \
    -o ${prefix}.gatk.bam

  ${params.samtools} index ${prefix}.gatk.bam
  """
}

process bedtools_count {

  storeDir "${baseDir}/BedtoolsCount/${prefix}"
  scratch true
  stageInMode 'copy'
  executor 'sge'
  clusterOptions '-l h_vmem=10G -pe smp 1 -l h_rt=12:00:00 -l os=rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_bedtools_count

  output:
    file('*.genome.counts.txt') into bedtools_count_out
    file('*.Trp53.counts.txt') into bedtools_count_trp53_out
    val prefix into samples_combine_bedtools

  """
  ~/chm2059/lib/bedtools2/bin/bedtools coverage \
    -a ${params.bedtools_intervals} \
    -b ${bam_file} \
    -counts -sorted \
    -g ${params.bedtools_genome} > ${prefix}.genome.counts.txt

  ~/chm2059/lib/bedtools2/bin/bedtools coverage \
    -a ${params.bedtools_intervals_trp53} \
    -b ${bam_file} \
    -counts -sorted \
    -g ${params.bedtools_genome} > ${prefix}.Trp53.counts.txt
  """
}

process combine_bedtools_count {

  storeDir "${baseDir}/BedtoolsCount/"

  input:
    file results from bedtools_count_out.toList()
    file results_trp53 from bedtools_count_trp53_out.toList()
    val all_samples from samples_combine_bedtools.toList()

  output:
    file('*.counts.csv') into combine_bedtools_count_out

  """
  python ${params.seqpy}/aggregate_bedtoolsCount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out BedtoolsCount.counts.csv

  python ${params.seqpy}/aggregate_bedtoolsCount.py \
    --files ${results_trp53} \
    --samples ${all_samples.join(" ")} \
    --out BedtoolsCount.Trp53.counts.csv
  """
}

process CollectInsertSizeMetrics {
  storeDir "${baseDir}/CollectInsertSizeMetrics/${prefix}"
  scratch true
  stageInMode 'copy'
  executor 'sge'
  maxForks 10
  clusterOptions '-l h_vmem=24G -pe smp 1 -l h_rt=24:00:00 -l os=rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_insert_size

  output:
    set prefix, file('*CollectInsertSizeMetrics.txt'), file('*pdf') into CollectInsertSizeMetrics_out

  """
  ${params.java} -Xmx16g -jar \
    ${params.picard} CollectInsertSizeMetrics \
    VALIDATION_STRINGENCY=LENIENT \
    H=${prefix}.pdf \
    I=${bam_file} \
    O=${prefix}.CollectInsertSizeMetrics.txt
  """
}

process dexseq_count {

  storeDir "${baseDir}/DEXSeqCount/${prefix}"
  scratch true
  stageInMode 'copy'
  executor 'sge'
  maxForks 10
  clusterOptions '-l h_vmem=2G -pe smp 4 -l h_rt=24:00:00 -l os=rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_dexseq_count

  output:
    file('*.dexseq.txt') into dexseq_count_out
    val prefix into samples_combine_dexseq

  """
  ${params.samtools} view -q 30 ${bam_file} | \
    ~/chm2059/lib/getnh | \
    python ~/chm2059/lib/DEXSeq/inst/python_scripts/dexseq_count.py \
    -p yes \
    -s no \
    -f sam \
    -r pos \
    ~/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.dexseq.gff \
    - \
    ${prefix}.dexseq.txt
  """
}

process combine_dexseq_count {

  storeDir "${baseDir}/DEXSeqCount/"

  input:
    file results from dexseq_count_out.toList()
    val all_samples from samples_combine_dexseq.toList()

  output:
    file('DEXseq.counts.csv') into combine_dexseq_count_out

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out DEXseq.counts.csv
  """
}

process gatk_coverage {

  storeDir "${baseDir}/GATKCoverageQ30/${prefix}"
  scratch true
  executor 'sge'
  stageInMode 'copy'
  maxForks 10
  clusterOptions '-l h_vmem=4G -pe smp 4 -l h_rt=24:00:00 -l os=rhel5.4|rhel6.3'

  input:
    set prefix, file(bam), file(bam_bai) from gatk_out_gatk_coverage

  output:
    val prefix into gatk_coverage_out_sample_1
    val prefix into gatk_coverage_out_sample_2
    file('*_summary') into gatk_coverage_out_summary
    file('*_statistics') into gatk_coverage_out_statistics
    file('*_interval_summary') into gatk_coverage_out_interval_summary
    file('*_interval_statistics') into gatk_coverage_out_interval_statistics
    file('*_cumulative_coverage_counts') into gatk_coverage_out_cumulative_coverage_counts
    file('*_cumulative_coverage_proportions') into gatk_coverage_out_cumulative_coverage_proportions

    script:
      if (prefix in initial_brca1_flfl)
        """
        ${params.samtools} view -h -q 30 \
          ${bam} \
          | awk \' \$7==\"=\" || \$1 ~ /^@/ { print \$0 }\' \
          | ${params.samtools} view -h -bS - > tmp.bam
        ${params.samtools} index tmp.bam

        ${params.java} -jar ${params.gatk} \
          -T DepthOfCoverage \
          -R ${params.reference} \
          -o ${prefix} \
          -L ${params.capture_file} \
          -ct 10 -ct 15 -ct 20 -ct 30 -ct 40 -ct 50 \
          -I tmp.bam
        """
      else
        """
        ${params.samtools} view -h -q 30 \
          ${bam} \
          | awk \' \$7==\"=\" || \$1 ~ /^@/ { print \$0 }\' \
          | ${params.samtools} view -h -bS - > tmp.bam
        ${params.samtools} index tmp.bam

        ${params.java} -jar ${params.gatk} \
          -T DepthOfCoverage \
          -R ${params.reference} \
          -o ${prefix} \
          -L ${params.capture_file_sureselect} \
          -ct 10 -ct 15 -ct 20 -ct 30 -ct 40 -ct 50 \
          -I tmp.bam
        """
}

gatk_coverage_combine_sureselect = Channel.create()
gatk_coverage_combine_nimblegen = Channel.create()

gatk_coverage_out_sample_1
  .map{
    "${baseDir}/GATKCoverageQ30/" + it
  }
  .merge(gatk_coverage_out_sample_2) {t1,t2 -> [t1,t2]}
  .choice(gatk_coverage_combine_nimblegen,gatk_coverage_combine_sureselect) { a -> a[1] in initial_brca1_flfl ? 0 : 1 }

gatk_coverage_combine_sureselec_sample = Channel.create()
gatk_coverage_combine_sureselect_dir = Channel.create()

gatk_coverage_combine_nimblegen_sample = Channel.create()
gatk_coverage_combine_nimblegen_dir = Channel.create()

gatk_coverage_combine_sureselect
  .separate(gatk_coverage_combine_sureselect_dir,gatk_coverage_combine_sureselec_sample) {a -> [a[0],a[1]]}

gatk_coverage_combine_nimblegen
  .separate(gatk_coverage_combine_nimblegen_dir,gatk_coverage_combine_nimblegen_sample) {a -> [a[0],a[1]]}


process gatk_coverage_combine_sureselect {

  storeDir "${baseDir}/GATKCoverageQ30/"

  input:
    val all_dirs from gatk_coverage_combine_sureselect_dir.toList()
    val all_samples from gatk_coverage_combine_sureselec_sample.toList()

  output:
    file('GATKCoverage.sureselect*') into gatk_coverage_combine_out_sureselect

  """
  python ${params.seqpy}/aggregate_gatk_coverage.py \
    --sample_dir ${all_dirs.join(" ")} \
    --samples ${all_samples.join(" ")} \
    --prefix GATKCoverage.sureselect
  """
}

process gatk_coverage_combine_nimblegen {

  storeDir "${baseDir}/GATKCoverageQ30/"

  input:
    val all_dirs from gatk_coverage_combine_nimblegen_dir.toList()
    val all_samples from gatk_coverage_combine_nimblegen_sample.toList()

  output:
    file('GATKCoverage.nimblegen*') into gatk_coverage_combine_out_nimblegen

  """
  python ${params.seqpy}/aggregate_gatk_coverage.py \
    --sample_dir ${all_dirs.join(" ")} \
    --samples ${all_samples.join(" ")} \
    --prefix GATKCoverage.nimblegen
  """
}

process samtools_flagstat_gatk {

  maxForks 10
  storeDir "${baseDir}/SamtoolsFlagstatGATK/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_flagstat

  output:
    file('*flagstat') into flagstat_out_gatk

  """
  ${params.samtools} flagstat ${bam_file} > ${prefix}.flagstat
  """
}

flagstats_gatk = flagstat_out_gatk.toList()

process combine_flagstat_gatk {

  storeDir "${baseDir}/SamtoolsFlagstatGATK"

  input:
    file flagstats_gatk
  output:
    file('flagstats.csv') into flagstat_gatk_combine_out

  """
  python ${params.flagstat} \
    --flagstat ${flagstats_gatk} \
    --out flagstats.csv
  """
}

process varscan_snps {

  storeDir "${baseDir}/Varscan/${prefix}"
  scratch true
  executor 'sge'
  maxForks 10
  stageInMode 'copy'
  clusterOptions '-l h_vmem=12G -pe smp 2 -l h_rt=48:00:00 -l os=rhel5.4|rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_varscan_snp

  output:
    file('*snps.vcf') into varscan_snp_out

  """
  rsync -L ${params.reference} ref.fa
  ${params.samtools} mpileup \
    -d10000 -f ref.fa \
    -q 15 \
    tmp.bam \
    | ${params.java} -Xmx16g -jar ${params.varscan} mpileup2snp \
    --p-value 0.05 \
    --output-vcf > \
    ${prefix}.varscan.snps.vcf
  """
}

process varscan_indel {

  storeDir "${baseDir}/Varscan/${prefix}"
  scratch true
  executor 'sge'
  maxForks 10
  stageInMode 'copy'
  clusterOptions '-l h_vmem=12G -pe smp 2 -l h_rt=48:00:00 -l os=rhel5.4|rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_varscan_indel

  output:
    file('*indels.vcf') into varscan_indel_out

  """
  rsync -L ${params.reference} ref.fa
  ${params.samtools} mpileup \
    -d10000 -f ref.fa \
    -q 15 \
    tmp.bam \
    | ${params.java} -Xmx16g -jar ${params.varscan} mpileup2indel \
    --p-value 0.05 \
    --output-vcf > \
    ${prefix}.varscan.indels.vcf
  """
}

def getFastqFiles( Path directory ) {

  def prefix = directory.name.split('_')[0]

  R1_files = []
  R2_files = []
  directory.eachFile{
    if (it.name.contains('R1')) {
      R1_files << it
    } else if (it.name.contains('R2')){
      R2_files << it
    }
  }
  R1_files.sort()
  R2_files.sort()
  return [prefix, R1_files, R2_files]
}
