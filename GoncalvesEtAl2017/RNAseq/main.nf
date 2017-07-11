#!/usr/bin/env nextflow

basedir = "$baseDir"

params.reference = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.dbsnp = '/home/chm2059/chm2059/data/refdata/mm10/mgp.v4.snps.dbSNP.vcf'
params.gtf = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.gtf'
params.dexseqgff = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.dexseq.gff'

params.star = '/home/chm2059/chm2059/lib/STAR-STAR_2.4.1d/bin/Linux_x86_64_static/STAR'
params.starRef = '/home/chm2059/chm2059/data/refdata/GRCm38/star100'
params.outFilterMismatchNmax = 5
params.star_threads = 6

params.dexseq = '/home/chm2059/chm2059/lib/DEXSeq/python_scripts/dexseq_count.py'
params.fastqscreen = '/home/chm2059/chm2059/lib/fastq_screen_v0.10.0/fastq_screen'
params.fastqscreen_conf = '/home/chm2059/chm2059/lib/fastq_screen_v0.10.0/fastq_screen.conf'
params.java = '/home/chm2059/chm2059/lib/jre1.8.0_25/bin/java'
params.fastqc = '/home/chm2059/chm2059/lib/FastQC/fastqc'
params.seqpy = '/home/chm2059/chm2059/lib/seqpy/bin/'
params.samtools = '/home/chm2059/chm2059/lib/samtools-1.1/samtools'
params.picard = '/home/chm2059/chm2059/lib/picard-tools-1.137/picard.jar'
params.htseq = '~/.local/bin/htseq-count'
params.flagstat = '~/chm2059/lib/seqpy/bin/flagstat.py'
params.cufflinks = '/home/chm2059/chm2059/lib/cufflinks-2.2.1.Linux_x86_64/cufflinks'

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()
fqfiles4 = Channel.create()

Channel
  .fromFilePairs("/home/chm2059/chm2059/0117_marcus/data/Marcus4376_2017_01_10/*_L00*_R1_001.fastq.gz", size: 1)
  .map{[it[0].split('_')[0],it[1][0]]}
  .separate(fqfiles, fqfiles2, fqfiles3, fqfiles4){it->[it,it,it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 6

  input:
    set prefix, file(read1) from fqfiles

  output:
    set file('*zip') into fastqc_results
    set prefix into samples

  """
  zcat ${read1} | ${params.fastqc} stdin --outdir=.
  mv stdin_fastqc.zip ${prefix}.zip
  """
}

process fastqc_screen {

  storeDir "${baseDir}/fastqc_screen/${prefix}"
  executor 'sge'
  clusterOptions '-l h_vmem=10G -pe smp 1 -l os=rhel6.3'

  input:
    set prefix, file(read1) from fqfiles4

  output:
    set file('*html'), file('*screen.txt') into fastqc_screen_results

  """
  ${params.fastqscreen} --aligner bwa \
    --conf ${params.fastqscreen_conf} \
    ${read1}
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
  clusterOptions '-l h_vmem=10G -pe smp 6 -l os=rhel6.3'

  storeDir "${baseDir}/STAR/${prefix}"

  input:
    set prefix, file(read1) from fqfiles3

  output:
    set prefix, file('*bam'), file('*bam.bai') into star_out
    set prefix, file('*bam'), file('*bam.bai') into star_out2
    set prefix, file('*bam'), file('*bam.bai') into star_out3
    set prefix, file('*bam'), file('*bam.bai') into star_out4
    set file('*Log.final.out') into star_log
    set file('*Log.progress.out') into star_log_process
    set file('*Log.out') into star_log_out
    set file('*Chimeric.out.junction') into star_chim_junction_out
    set file('*Chimeric.out.sam') into star_chim_sam_out
    set file('*bam') into bam_junction_count
    set file('*bam.bai') into bam_bai_junction_count
    set prefix into samples_star_qc
    set prefix into samples_junction_count


  script:
  read1=read1.toString().replaceAll(/ /,",")

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
    --readFilesIn ${read1} \
    --outFileNamePrefix ./1pass/tmp.

  ${params.star} \
    --runMode genomeGenerate \
    --genomeDir ./star_2pass \
    --genomeFastaFiles ${params.reference} \
    --sjdbFileChrStartEnd ./1pass/tmp.SJ.out.tab \
    --sjdbOverhang 100 \
    --runThreadN ${params.star_threads}

  ${params.star} \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN ${params.star_threads} \
    --outSAMstrandField intronMotif \
    --genomeDir ./star_2pass \
    --readFilesCommand zcat \
    --readFilesIn ${read1} \
    --outFileNamePrefix ${prefix}. \
    --chimOutType SeparateSAMold \
    --chimSegmentMin 1

   ${params.samtools} index ${prefix}*.bam
   rm -rf 1pass
   rm -rf star_2pass
 """
}

process samtools_flagstat {

  maxForks 8
  storeDir "${baseDir}/SamtoolsFlagstat/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out

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
  clusterOptions '-l h_vmem=3G -pe smp 6 -l os=rhel6.3'

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
    -q -p 6 \
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

  storeDir "${baseDir}/CufflinksSTAR/"

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
  executor 'sge'
  clusterOptions '-l h_vmem=16G -pe smp 1 -l h_rt=24:00:00 -l os=rhel6.3'

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out3

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
    file('HTSeq.gene.counts.csv') into combine_htseq_out

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out HTSeq.gene.counts.csv
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
