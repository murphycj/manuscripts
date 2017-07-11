#!/usr/bin/env nextflow

PE75_samples = [
  'HL-167-2','HL115-2374','HL158','HL160','HL161','HL162','HL163','HL164',
  'HL165','HL166','HL168','HL169','HL170','HL171','HL172','HL173','HL174',
  'HL175','HL176','HL177','HL116-2453','HL117-2598'
]

params.fastqs = "/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/data/RNAseq/Sample_*"

basedir = "$baseDir"

params.reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.reference_dict = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.dict'
params.reference_fai = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'
params.dbsnp = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/mm10/mgp.v4.snps.dbSNP.vcf'
params.dexseq_gff = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.dexseq.gff'
params.gtf = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.gtf'

params.star = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/STAR-STAR_2.4.1d/bin/Linux_x86_64_static/STAR'
params.starRef = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/star'
params.outFilterMismatchNmax = 5
params.star_threads = 4

params.fusioncatcher = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/fusioncatcher/bin/fusioncatcher'
params.fusioncatcher_ref = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/mm10/fusioncatcher/'
params.fusioncatcher_threads = 4

params.varscan = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/VarScan.v2.3.9.jar'
params.java = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/jre1.8.0_25/bin/java'
params.fastqc = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/FastQC/fastqc'
params.seqpy = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/'
params.samtools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/samtools-1.3.1/samtools'
params.picard = '~/chm2059/lib/picard-tools-1.137/picard.jar'
params.gatk = '~/chm2059/lib/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
params.picard = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/picard-tools-1.137/picard.jar'
params.htseq = '~/.local/bin/htseq-count'
params.flagstat = '~/chm2059/lib/seqpy/bin/flagstat.py'
params.cufflinks = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/cufflinks-2.2.1.Linux_x86_64/cufflinks'

log.info 'RNA-seq pipeline'

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()

Channel
  .fromPath(params.fastqs, type: 'dir')
  .map { path ->
       (prefix, R1_files, R2_files, rl) = getFastqFiles(path)
       tuple(prefix, R1_files, R2_files, rl)
  }
  .separate(fqfiles, fqfiles2, fqfiles3){it->[it,it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 10

  input:
    set prefix, file(read1), file(read2), rl from fqfiles

  output:
    file('*zip') into fastqc_results
    val prefix into samples

  """
  zcat ${read1} ${read2} | ${params.fastqc} stdin --outdir=.
  mv stdin_fastqc.zip ${prefix}.zip
  """
}

process combine_fastqc_results {

  storeDir "${baseDir}/fastqc"

  input:
    file results from fastqc_results.toList()
    val all_samples from samples.toList()

  output:
    file('HL_RNAseq_Fastqc.csv') into fastqc_combine_output

  """
  python ${params.seqpy}/aggregate_fastqc.py \
    --out HL_RNAseq_Fastqc.csv \
    --files ${results} \
    --samples ${all_samples.join(" ")}
  """
}

process fusioncatcher {
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=6G -pe smp 4 -l athena=true'

  storeDir "${baseDir}/fusioncatcher/"

  input:
    set prefix, file(read1), file(read2), rl from fqfiles2

  output:
    file(prefix) into fusioncatcher_out

  """
  cat ${read1} > ${prefix}.R1.fastq.gz
  cat ${read2} > ${prefix}.R2.fastq.gz

  ${params.fusioncatcher} \
    -p ${params.fusioncatcher_threads} \
    -d ${params.fusioncatcher_ref} \
    -i ${prefix}.R1.fastq.gz,${prefix}.R2.fastq.gz \
    -o ${prefix}
  """
}

process process_gatk {
  executor 'sge'
  scratch true
  maxForks 20
  clusterOptions '-l h_vmem=10G -pe smp 4 -l h_rt=96:00:00 -l athena=true'

  storeDir "${baseDir}/STARGATK/${prefix}"

  input:
    set prefix, file(read1), file(read2), rl from fqfiles3

  output:
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out2
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_insert_size
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out3
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_varscan_snp
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_varscan_indel
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_mutect
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into star_out_htseq_exons
    set file('*Log.progress.out') into star_out_progress
    set file('*Log.final.out') into star_log
    set file('*SJ.out.tab') into star_out_SJ

  script:
  read1=read1.toString().replaceAll(/ /,",")
  read2=read2.toString().replaceAll(/ /,",")

  """
  mkdir 1pass
  mkdir star_2pass
  ${params.star} \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN ${params.star_threads} \
    --outSAMstrandField intronMotif \
    --outFilterMismatchNmax ${params.outFilterMismatchNmax} \
    --genomeDir ${params.starRef}${rl} \
    --readFilesCommand zcat \
    --readFilesIn ${read1} ${read2} \
    --outFileNamePrefix ./1pass/tmp.

  ${params.star} \
    --runMode genomeGenerate \
    --genomeDir ./star_2pass \
    --genomeFastaFiles ${params.reference} \
    --sjdbFileChrStartEnd ./1pass/tmp.SJ.out.tab \
    --sjdbOverhang ${rl} \
    --runThreadN ${params.star_threads}

  ${params.star} \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN ${params.star_threads} \
    --outSAMstrandField intronMotif \
    --genomeDir ./star_2pass \
    --readFilesCommand zcat \
    --readFilesIn ${read1} ${read2} \
    --outFileNamePrefix ${prefix}.

   ${params.samtools} index ${prefix}.Aligned.sortedByCoord.out.bam
   rm -rf 1pass
   rm -rf star_2pass

  ${params.java} -Xmx10g -jar \
    ${params.picard} AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=LENIENT \
    I=${prefix}.Aligned.sortedByCoord.out.bam \
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
    REMOVE_DUPLICATES=false \
    M=duplicates.bam

  ${params.samtools} index tmp2.bam
  rm tmp.bam*

  ${params.java} -jar \
    ${params.gatk} \
    -T SplitNCigarReads \
    -R ${params.reference} \
    -I tmp2.bam \
    -o tmp.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

  ${params.samtools} index tmp.bam
  rm tmp2.bam*

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T RealignerTargetCreator \
    -R ${params.reference} \
    -I tmp.bam \
    -U ALLOW_N_CIGAR_READS \
    -o forIndelRealigner.intervals \
    -nt 4

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T IndelRealigner \
    -R ${params.reference} \
    -I tmp.bam \
    -targetIntervals forIndelRealigner.intervals \
    -o tmp2.bam

  ${params.samtools} index tmp2.bam

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T BaseRecalibrator \
    -R ${params.reference} \
    -I tmp2.bam \
    -o recal_data.table \
    -knownSites ${params.dbsnp} \
    -nct 4

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T PrintReads \
    -R ${params.reference} \
    -I tmp2.bam \
    --BQSR recal_data.table \
    -o tmp3.bam

  ${params.samtools} index tmp3.bam

  python ${params.seqpy}/filter_bad_cigar.py \
    --infile tmp3.bam \
    --out ${prefix}.gatk.bam

  ${params.samtools} index ${prefix}.gatk.bam
  """
}

process CollectInsertSizeMetrics {
  storeDir "${baseDir}/CollectInsertSizeMetrics/${prefix}"
  scratch true
  stageInMode 'copy'
  executor 'sge'
  maxForks 10
  clusterOptions '-l h_vmem=24G -pe smp 1 -l h_rt=24:00:00 -l athena=true'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_insert_size

  output:
    set prefix, file('*CollectInsertSizeMetrics.txt') into CollectInsertSizeMetrics_out

  """
  ${params.java} -Xmx16g -jar \
    ${params.picard} CollectInsertSizeMetrics \
    VALIDATION_STRINGENCY=LENIENT \
    H=${prefix}.pdf \
    I=${bam_file} \
    O=${prefix}.CollectInsertSizeMetrics.txt
  """
}

process samtools_flagstatGATK {

  maxForks 8
  storeDir "${baseDir}/SamtoolsFlagstatGATK/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out3

  output:
    file('*flagstat') into flagstat_outGATK

  """
  ${params.samtools} flagstat ${bam_file} > ${prefix}.flagstat

  """
}

process combine_flagstatGATK {

  storeDir "${baseDir}/SamtoolsFlagstatGATK"

  input:
    file(flagstat_files) from flagstat_outGATK.toList()
  output:
    file('flagstats.csv') into flagstat_combine_out_GATK

  """
  python ${params.seqpy}/aggregate_flagstat.py \
    --flagstat ${flagstat_files} \
    --out flagstats.csv
  """
}

process cufflinks_gatk {
  scratch true
  maxForks 10

  storeDir "${baseDir}/Cufflinks/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out
  output:
    set prefix, file('*genes.fpkm_tracking'), file('*isoforms.fpkm_tracking'), file('*skipped.gtf'), file('*transcripts.gtf') into cufflinks_out_gatk
    file('*genes.fpkm_tracking') into cufflinks_out_genes_gatk
    file('*isoforms.fpkm_tracking') into cufflinks_out_isoforms_gatk
    val prefix into samples_cufflinks_gatk

  """
  ${params.cufflinks} \
    -q -p 4 \
    -o . \
    -b ${params.reference} \
    -G ${params.gtf} \
    ${bam_file}
  mv genes.fpkm_tracking ${prefix}.genes.fpkm_tracking
  mv isoforms.fpkm_tracking ${prefix}.isoforms.fpkm_tracking
  mv skipped.gtf ${prefix}.skipped.gtf
  mv transcripts.gtf ${prefix}.transcripts.gtf
  """
}

process combine_cufflinks {

  storeDir "${baseDir}/Cufflinks/"

  input:
    file genes from cufflinks_out_genes_gatk.toList()
    file isoforms from cufflinks_out_isoforms_gatk.toList()
    val all_samples from samples_cufflinks_gatk.toList()

  output:
    file('isoforms-fpkm.csv') into combine_cufflinks_isoforms_gatk
    file('genes-fpkm.csv') into combine_cufflinks_genes_gatk

  """
  python ${params.seqpy}/aggregate_fpkm.py \
    --files ${genes} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out genes-fpkm.csv

  python ${params.seqpy}/aggregate_fpkm.py \
    --files ${isoforms} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out isoforms-fpkm.csv

  """
}

process convert_expression_files {

  storeDir "${baseDir}/Cufflinks/"

  input:
    set genes from combine_cufflinks_genes_gatk

  output:
    file('genes-symbols-fpkm.csv') into genes_symbols_gatk
    file('genes-humanSymbols-fpkm.csv') into genes_symbols_human_gatk

  """
  Rscript ${params.seqpy}/ensembl2symbol_mouse.R \
    -fpkm ${genes} \
    -out genes-symbols-fpkm.csv

  python ${params.seqpy}/mouseSymbol2human.py \
    --infile genes-symbols-fpkm.csv \
    --out genes-humanSymbols-fpkm.csv
  """
}


process htseq_readsGATK {

  storeDir "${baseDir}/HTSeqCount/${prefix}"
  maxForks 30

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out2

  output:
    file('*count') into htseq_reads_outGATK
    val prefix into samples_combine_htseqGATK

  """
  python ${params.htseq} \
    -s no \
    -f bam \
    ${bam_file} \
    ${params.gtf} > ${prefix}.count
  """
}

process combine_htseqGATK {

  storeDir "${baseDir}/HTSeqCount/"

  input:
    file results from htseq_reads_outGATK.toList()
    val all_samples from samples_combine_htseqGATK.toList()

  output:
    file('HTSeq.gene.counts.csv') into combine_htseq_outGATK

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out HTSeq.gene.counts.csv
  """
}


process htseq_reads_exons {

  storeDir "${baseDir}/HTSeqCount/${prefix}"
  maxForks 35

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out_htseq_exons

  output:
    file('*.exoncount.txt') into htseq_reads_exons_out
    val prefix into samples_combine_htseq_exons

  """
  python ${params.htseq} \
    -s no \
    -f bam \
    -t exon \
    -i exon_id \
    ${bam_file} \
    ${params.gtf} > ${prefix}.exoncount.txt
  """
}

process combine_htseq_exons {

  storeDir "${baseDir}/HTSeqCount/"

  input:
    file results from htseq_reads_exons_out.toList()
    val all_samples from samples_combine_htseq_exons.toList()

  output:
    file('HTSeq.exon.counts.csv') into combine_htseq_exons_out

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out HTSeq.exon.counts.csv
  """
}

process varscan_snps {

  storeDir "${baseDir}/Varscan/${prefix}"
  scratch true
  executor 'sge'
  stageInMode 'copy'
  maxForks 10
  clusterOptions '-l h_vmem=12G -pe smp 2 -l h_rt=48:00:00 -l athena=true'

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_varscan_snp

  output:
    file('*.snps.vcf') into varscan_snp_out

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
  clusterOptions '-l h_vmem=12G -pe smp 2 -l h_rt=48:00:00 -l athena=true'

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

  def prefix = directory.name.split('_')[1]

  def rl = '50'
  if (prefix in PE75_samples) {
    rl = '74'
  } else {
    rl = '50'
  }

  R1_files = []
  R2_files = []
  directory.eachFile{
    if (it.name.contains('R1')) {
      R1_files << it
    } else if (it.name.contains('R2')) {
      R2_files << it
    }
  }
  R1_files.sort()
  R2_files.sort()
  return [prefix, R1_files, R2_files, rl]
}
