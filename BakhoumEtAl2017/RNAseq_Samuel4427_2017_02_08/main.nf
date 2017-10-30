#!/usr/bin/env nextflow


params.fastqs = "/home/chm2059/chm2059/0816_sam/data/Samuel4427_2017_02_08/*"

basedir = "$baseDir"

params.reference = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.reference_dict = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.dict'
params.reference_fai = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'
params.dbsnp = '/home/chm2059/chm2059/data/refdata/GRCh38/common_all_20160407.vcf'
params.gtf = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.88.chr.gtf'

params.starRef = '/home/chm2059/chm2059/data/refdata/GRCh38/star100'
params.starOverhang = 100
params.outFilterMismatchNmax = 5
params.star_threads = 12

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()
fqfiles4 = Channel.create()

fastqs = Channel
  .fromFilePairs("/athena/elementolab/scratch/chm2059/from_dat02/chm2059/0816_sam/data/xenome/*/*graft.fastq.gz", size: 1)
  .map{[it[0].split('_')[0],it[1][0]]}

Channel
  .fromFilePairs("/home/chm2059/chm2059/0816_sam/data/Samuel4427_2017_02_08/*_L00*_R1_001.fastq.gz", size: 1)
  .map{[it[0].split('_')[0],it[1][0]]}
  .filter{(it[0] in['11','12','13','14','15','16','17','18'])}
  .concat(fastqs)
  .separate(fqfiles, fqfiles2, fqfiles3, fqfiles4){it->[it,it,it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 6

  input:
    set prefix, file(read1) from fqfiles

  output:
    file('*zip') into fastqc_results
    val prefix into samples

  """
  zcat ${read1} | ${params.fastqc} stdin --outdir=.
  mv stdin_fastqc.zip ${prefix}.zip
  """
}

process combine_fastqc_results {

  storeDir "${baseDir}/fastqc"

  input:
    file results from fastqc_results.toList()
    val all_samples from samples.toList()

  output:
    file('RNAseq_Fastqc.csv') into fastqc_combine_output

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
  clusterOptions '-l h_vmem=4G -pe smp 12 -l h_rt=96:00:00 -l athena=true'

  storeDir "${baseDir}/STAR/${prefix}"

  input:
    set prefix, file(read1) from fqfiles3

  output:
    set prefix, file('*bam'), file('*bam.bai') into star_out
    set prefix, file('*bam'), file('*bam.bai') into star_out2
    set prefix, file('*bam'), file('*bam.bai') into star_out3
    set prefix, file('*bam'), file('*bam.bai') into star_out4
    file('*Log.final.out') into star_log
    file('*Log.progress.out') into star_log_process
    file('*Log.out') into star_log_out
    file('*Chimeric.out.junction') into star_chim_junction_out
    file('*Chimeric.out.sam') into star_chim_sam_out
    file('*bam') into bam_junction_count
    file('*bam.bai') into bam_bai_junction_count
    val prefix into samples_star_qc
    val prefix into samples_junction_count
    file('.command*') into star_scripts


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

  maxForks 12
  storeDir "${baseDir}/SamtoolsFlagstat/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file) from star_out3

  output:
    file('*flagstat') into flagstat_out

  """
  ${params.samtools} flagstat ${bam_file} > ${prefix}.flagstat

  """
}

process combine_flagstat {

  storeDir "${baseDir}/SamtoolsFlagstat"

  input:
    file(flagstat_files) from flagstat_out.toList()
  output:
    file('flagstats.csv') into flagstat_combine_out

  """
  python ${params.seqpy}/aggregate_flagstat.py \
    --flagstat ${flagstat_files} \
    --out flagstats.csv
  """
}

process cufflinks_star {
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=2G -pe smp 12 -l h_rt=96:00:00 -l athena=true'

  storeDir "${baseDir}/CufflinksSTAR/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file) from star_out2

  output:
    set prefix, file('*genes.fpkm_tracking'), file('*isoforms.fpkm_tracking'), file('*skipped.gtf'), file('*transcripts.gtf') into cufflinks_out
    file('*genes.fpkm_tracking') into cufflinks_out_genes
    val prefix into samples_cufflinks

  """
  ${params.cufflinks} \
    -q -p 12 \
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
    val all_samples from samples_cufflinks.toList()

  output:
    file('genes-fpkm.csv') into combine_cufflinks_genes

  """
  python /home/chm2059/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${genes} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out genes-fpkm.csv
  """
}

process htseq_reads {

  storeDir "${baseDir}/HTSeqCount/${prefix}"
  maxForks 20

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out

  output:
    file('*count') into htseq_reads_out
    val prefix into samples_combine_htseq

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
