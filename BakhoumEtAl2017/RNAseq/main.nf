#!/usr/bin/env nextflow


params.fastqs = "/home/chm2059/chm2059/0816_sam/data/RNAseq/*"

basedir = "$baseDir"

params.reference = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.reference_dict = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.dict'
params.reference_fai = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
params.dbsnp = '/home/chm2059/chm2059/data/refdata/GRCh38/common_all_20160407.vcf'
params.gtf = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.79.gtf'

params.star = '/home/chm2059/chm2059/lib/STAR-STAR_2.4.1d/bin/Linux_x86_64_static/STAR'
params.starRef = '/home/chm2059/chm2059/data/refdata/GRCh38/star'
params.starOverhang = 74
params.outFilterMismatchNmax = 5
params.star_threads = 4

params.java = '/home/chm2059/chm2059/lib/jre1.8.0_25/bin/java'

params.varscan = '/home/chm2059/chm2059/lib/VarScan.v2.3.9.jar'

params.fastqc = '/home/chm2059/chm2059/lib/FastQC/fastqc'

params.seqpy = '/home/chm2059/chm2059/lib/seqpy/bin/'

params.samtools = '/home/chm2059/chm2059/lib/samtools-1.1/samtools'

params.picard = '~/chm2059/lib/picard-tools-1.137/picard.jar'

params.gatk = '~/chm2059/lib/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'

params.picard = '/home/chm2059/chm2059/lib/picard-tools-1.137/picard.jar'

params.htseq = '~/.local/bin/htseq-count'

params.flagstat = '~/chm2059/lib/seqpy/bin/flagstat.py'

params.cufflinks = '/home/chm2059/chm2059/lib/cufflinks-2.2.1.Linux_x86_64/cufflinks'

log.info 'RNA-seq pipeline'

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()

Channel
  .fromPath(params.fastqs, type: 'dir')
  .map { path ->
       (prefix, R1_files, R2_files) = getFastqFiles(path)
       tuple(prefix, R1_files, R2_files)
  }
  .filter{
    it[0]!='Summary'
  }
  .separate(fqfiles, fqfiles2, fqfiles3){it->[it,it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 6

  input:
    set prefix, file(read1), file(read2) from fqfiles

  output:
    set file('*zip') into fastqc_results
    set prefix into samples

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
    set file('RNAseq_Fastqc.csv') into fastqc_combine_output

  """
  python ${params.seqpy}/aggregate_fastqc.py \
    --out RNAseq_Fastqc.csv \
    --files ${results} \
    --samples ${all_samples.join(" ")}
  """
}

process star {
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=10G -pe smp 4'

  storeDir "${baseDir}/STAR/${prefix}"

  input:
    set prefix, file(read1), file(read2) from fqfiles3

  output:
    set prefix, file('*bam'), file('*bam.bai') into star_out
    set prefix, file('*bam'), file('*bam.bai') into star_out2
    set prefix, file('*bam'), file('*bam.bai') into star_out3
    set prefix, file('*bam'), file('*bam.bai') into star_out_varscan_snp
    set prefix, file('*bam'), file('*bam.bai') into star_out_varscan_indel
    set prefix, file('*bam'), file('*bam.bai') into star_out_gatk
    set file('*Log.final.out') into star_log
    set prefix into samples_star_qc
    set file('*SJ.out.tab') into star_out_SJ
    set file('*Log.progress.out') into star_out_progress


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
    --genomeDir ${params.starRef} \
    --readFilesCommand zcat \
    --readFilesIn ${read1} ${read2} \
    --outFileNamePrefix ./1pass/tmp.

  ${params.star} \
    --runMode genomeGenerate \
    --genomeDir ./star_2pass \
    --genomeFastaFiles ${params.reference} \
    --sjdbFileChrStartEnd ./1pass/tmp.SJ.out.tab \
    --sjdbOverhang ${params.starOverhang} \
    --runThreadN ${params.star_threads}

  ${params.star} \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN ${params.star_threads} \
    --outSAMstrandField intronMotif \
    --genomeDir ./star_2pass \
    --readFilesCommand zcat \
    --readFilesIn ${read1} ${read2} \
    --outFileNamePrefix ${prefix}.

   ${params.samtools} index ${prefix}*.bam
   rm -rf 1pass
   rm -rf star_2pass
 """
}

process samtools_flagstat {

  maxForks 8
  storeDir "${baseDir}/SamtoolsFlagstat/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file) from star_out3

  output:
    set file('*flagstat') into flagstat_out

  """
  ${params.samtools} flagstat ${bam_file} > ${prefix}.flagstat

  """
}

process combine_flagstat {

  storeDir "${baseDir}/SamtoolsFlagstat"

  input:
    set file(flagstat_files) from flagstat_out.toList()
  output:
    set file('flagstats.csv') into flagstat_combine_out

  """
  python ${params.flagstat} \
    --flagstat ${flagstat_files} \
    --out flagstats.csv
  """
}

process cufflinks_star {
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=6G -pe smp 2'

  storeDir "${baseDir}/CufflinksSTAR/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file) from star_out2
  output:
    set prefix, file('*genes.fpkm_tracking'), file('*isoforms.fpkm_tracking'), file('*skipped.gtf'), file('*transcripts.gtf') into cufflinks_out
    set file('*genes.fpkm_tracking') into cufflinks_out_genes
    set file('*isoforms.fpkm_tracking') into cufflinks_out_isoforms
    set prefix into samples_cufflinks

  """
  ${params.cufflinks} \
    -q -p 2 \
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
    file genes from cufflinks_out_genes.toList()
    file isoforms from cufflinks_out_isoforms.toList()
    val all_samples from samples_cufflinks.toList()

  output:
    set file('isoforms-fpkm.csv') into combine_cufflinks_isoforms
    set file('genes-fpkm.csv') into combine_cufflinks_genes

  """
  python /home/chm2059/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${genes} \
    --samples ${all_samples.join(" ")} \
    --remove_duplicates \
    --out genes-fpkm.csv

  python /home/chm2059/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${isoforms} \
    --samples ${all_samples.join(" ")} \
    --remove_duplicates \
    --out isoforms-fpkm.csv
  """
}

process htseq_reads {

  storeDir "${baseDir}/HTSeqCount/${prefix}"
  maxForks 8

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out

  output:
    set file('*count') into htseq_reads_out
    set prefix into samples_combine_htseq

  """
  python ${params.htseq} \
    -s no \
    -f bam \
    ${bam_file} \
    ${params.gtf} > ${prefix}.count
  """
}

process combine_htseq {

  storeDir "${baseDir}/HTSeqCount/"

  input:
    file results from htseq_reads_out.toList()
    val all_samples from samples_combine_htseq.toList()

  output:
    file('HTSeq.gene.counts.csv') into combine_htseq_out_count
    file('HTSeq.geneSymbols.counts.csv') into combine_htseq_out_countSymbol
    file('HTSeq.human.geneSymbols.counts.csv') into combine_htseq_out_countHumanSymbol

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out HTSeq.gene.counts.csv

  Rscript ${params.seqpy}/ensembl2symbol_mouse.R \
    -data HTSeq.gene.counts.csv \
    -out HTSeq.geneSymbols.counts.csv

  python ${params.seqpy}/mouseSymbol2human.py \
    --infile HTSeq.geneSymbols.counts.csv \
    --out HTSeq.human.geneSymbols.counts.csv
  """
}

process combine_star_qc {

  storeDir "${baseDir}/STAR/"

  input:
    file log_files from star_log.toList()
    val all_samples from samples_star_qc.toList()

  output:
    file('STAR.QC.csv') into star_qc_out

  """
  python ${params.seqpy}/aggregate_star_QC.py \
    --files ${log_files} \
    --samples ${all_samples.join(" ")} \
    --out STAR.QC.csv
  """
}


process process_gatk {
  executor 'sge'
  scratch true
  clusterOptions '-l h_vmem=6G -pe smp 4 -l h_rt=72:00:00 -l os=rhel5.4|rhel6.3'

  storeDir "${baseDir}/STARGATK/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file) from star_out_gatk

  output:
    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out_mutect

  """
  ${params.java} -Xmx10g -jar \
    ${params.picard} AddOrReplaceReadGroups \
    VALIDATION_STRINGENCY=LENIENT \
    I=${bam_file} \
    O=tmp.bam \
    ID=test \
    LB=1 \
    PL=illumina \
    PU=1 \
    RGSM=test

  ${params.samtools} index tmp.bam

  ${params.java} -jar \
    ${params.gatk} \
    -T SplitNCigarReads \
    -R ${params.reference} \
    -I tmp.bam \
    -o split.bam \
    -rf ReassignOneMappingQuality \
    -RMQF 255 \
    -RMQT 60 \
    -U ALLOW_N_CIGAR_READS

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T BaseRecalibrator \
    -R ${params.reference} \
    -I split.bam \
    -o recal_data.table \
    -knownSites ${params.dbsnp} \
    -nct 4

  ${params.java} -Xmx10g -jar \
    ${params.gatk} \
    -T PrintReads \
    -R ${params.reference} \
    -I split.bam \
    --BQSR recal_data.table \
    -o ${prefix}.gatk.bam

  ${params.samtools} index ${prefix}.gatk.bam
  """
}

process mutect {

  storeDir "${baseDir}/Mutect/${prefix}"
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=12G -pe smp 6 -l h_rt=72:00:00 -l os=rhel5.4|rhel6.3'
  maxForks 30

  input:
    set prefix, file(bam_file), file(bam_index_file) from gatk_out_mutect

  output:
    set file('*vcf') into mutect_out

  """
  rsync -L ${params.reference} ref.fa
  rsync -L ${params.reference_dict} ref.dict
  rsync -L ${params.reference_fai} ref.fa.fai
  rsync -L ${bam_file} tmp.bam
  rsync -L ${bam_index_file} tmp.bam.bai
  ${params.java} -Xmx16g -jar \
    ${params.gatk} \
    -T MuTect2 \
    -R ref.fa \
    -I:tumor tmp.bam \
    -nct 6 \
    -o ${prefix}.vcf
  """
}

process varscan_snps {

  storeDir "${baseDir}/Varscan/${prefix}"
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=12G -pe smp 2  -l h_rt=72:00:00 -l os=rhel5.4|rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out_varscan_snp

  output:
    set file('*vcf') into varscan_snp_out

  """
  rsync -L ${params.reference} ref.fa
  rsync -L ${params.reference_dict} ref.dict
  rsync -L ${params.reference_fai} ref.fa.fai
  rsync -L ${bam_file} tmp.bam
  rsync -L ${bam_index_file} tmp.bam.bai
  ${params.samtools} mpileup \
    -d10000 -f ref.fa \
    tmp.bam \
    | ${params.java} -Xmx7g -jar ${params.varscan} mpileup2snp \
    --p-value 0.05 \
    --output-vcf > \
    ${prefix}.varscan.snps.vcf
  """
}

process varscan_indel {

  storeDir "${baseDir}/Varscan/${prefix}"
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=12G -pe smp 2  -l h_rt=72:00:00 -l os=rhel5.4|rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out_varscan_indel

  output:
    set file('*vcf') into varscan_indel_out

  """
  rsync -L ${params.reference} ref.fa
  rsync -L ${params.reference_dict} ref.dict
  rsync -L ${params.reference_fai} ref.fa.fai
  rsync -L ${bam_file} tmp.bam
  rsync -L ${bam_index_file} tmp.bam.bai
  ${params.samtools} mpileup \
    -d10000 -f ref.fa \
    tmp.bam \
    | ${params.java} -Xmx7g -jar ${params.varscan} mpileup2indel \
    --p-value 0.05 \
    --output-vcf > \
    ${prefix}.varscan.indels.vcf
  """
}

def getFastqFiles( Path directory ) {

  def prefix = directory.name

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
  return [prefix, R1_files, R2_files]
}
