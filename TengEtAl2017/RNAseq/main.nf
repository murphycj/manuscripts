#!/usr/bin/env nextflow

params.fastqs = "/home/chm2059/chm2059/dow_2016_5/data/Teng3670_2016_05_17/*"

basedir = "$baseDir"

params.reference = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.dbsnp = '/home/chm2059/chm2059/data/refdata/mm10/mgp.v4.snps.dbSNP.vcf'
params.gtf = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.gtf'
params.dexseqgff = '/home/chm2059/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.79.dexseq.gff'

params.dexseq = '/home/chm2059/chm2059/lib/DEXSeq/python_scripts/dexseq_count.py'

params.star = '/home/chm2059/chm2059/lib/STAR-STAR_2.4.1d/bin/Linux_x86_64_static/STAR'
params.starRef = '/home/chm2059/chm2059/data/refdata/GRCm38/star'
params.outFilterMismatchNmax = 5
params.star_threads = 4

params.fusioncatcher = '/home/chm2059/chm2059/lib/fusioncatcher/bin/fusioncatcher'
params.fusioncatcher_ref = '/home/chm2059/chm2059/data/refdata/mm10/fusioncatcher/'
params.fusioncatcher_threads = 4

params.java = '/home/chm2059/chm2059/lib/jre1.8.0_25/bin/java'
params.fastqc = '/home/chm2059/chm2059/lib/FastQC/fastqc'
params.seqpy = '/home/chm2059/chm2059/lib/seqpy/bin/'
params.samtools = '/home/chm2059/chm2059/lib/samtools-1.1/samtools'
params.picard = '~/chm2059/lib/picard-tools-1.137/picard.jar'
params.gatk = '~/chm2059/lib/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar'
params.picard = '/home/chm2059/chm2059/lib/picard-tools-1.137/picard.jar'
params.htseq = '~/.local/bin/htseq-count'
params.flagstat = '~/chm2059/lib/seqpy/bin/flagstat.py'
params.cufflinks = '/home/chm2059/chm2059/lib/cufflinks-2.2.1.Linux_x86_64/cufflinks'

fqfiles = Channel.create()
fqfiles2 = Channel.create()
fqfiles3 = Channel.create()

Channel
  .fromPath(params.fastqs, type: 'dir')
  .map { path ->
       (prefix, R1_files, rl) = getFastqFiles(path)
       tuple(prefix, R1_files, rl)
  }
  .filter{
    it[0]!='Summary'
  }
  .separate(fqfiles, fqfiles2, fqfiles3){it->[it,it,it]}


  deseq2_comparisons = Channel
    .from(
      [
        'APC_vs_WT',
        'APC',
        'WT',
        ['A1_plus','A2_plus','A5_plus','A6_plus'],
        ['A1','A2','A5','A6','NAIVE','LIPO']
      ],
      [
        'WT_vs_Rspondin',
        'Group2',
        'Group3',
        ['A1','A2','A5','A6','NAIVE','LIPO'],
        ['EN15d','EN4d','EN_M','EN_F','NAIVE_F']
      ],
      [
        'APC_vs_Rspondin',
        'APC',
        'Rspondin',
        ['A1_plus','A2_plus','A5_plus','A6_plus'],
        ['EN15d','EN4d','EN_M','EN_F','NAIVE_F']
      ]
    )

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 6

  input:
    set prefix, file(read1), rl from fqfiles

  output:
    set file('*zip') into fastqc_results
    set prefix into samples

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
    set prefix, file(read1), rl from fqfiles3

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
    --genomeDir ${params.starRef}${rl} \
    --readFilesCommand zcat \
    --readFilesIn ${read1} \
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
    --readFilesIn ${read1} \
    --outFileNamePrefix ${prefix}. \
    --chimOutType SeparateSAMold \
    --chimSegmentMin 1

   ${params.samtools} index ${prefix}*.bam
   rm -rf 1pass
   rm -rf star_2pass
 """
}

//process process_gatk {
//  executor 'sge'
//  scratch true
//  clusterOptions '-l h_vmem=6G -pe smp 4'

//  storeDir "${baseDir}/STARGATK/${prefix}"

//  input:
//    set prefix, file(bam_file), file(bam_bai_file) from star_out

//  output:
//    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out
//    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out2
//    set prefix, file('*.gatk.bam'), file('*.gatk.bam.bai') into gatk_out3

//  """
//  ${params.java} -Xmx10g -jar \
//    ${params.picard} AddOrReplaceReadGroups \
//    VALIDATION_STRINGENCY=LENIENT \
//    I=${bam_file} \
//    O=tmp.bam \
//    ID=test \
//    LB=1 \
//    PL=illumina \
//    PU=1 \
//    RGSM=test

//  ${params.samtools} index tmp.bam

//  ${params.java} -jar \
//    ${params.gatk} \
//    -T SplitNCigarReads \
//    -R ${params.reference} \
//    -I tmp.bam \
//    -o split.bam \
//    -rf ReassignOneMappingQuality \
//    -RMQF 255 \
//    -RMQT 60 \
//    -U ALLOW_N_CIGAR_READS

//  ${params.java} -Xmx10g -jar \
//    ${params.gatk} \
//    -T RealignerTargetCreator \
//    -R ${params.reference} \
//    -I split.bam \
//    -U ALLOW_N_CIGAR_READS \
//    -o forIndelRealigner.intervals \
//    -nt 4

//  ${params.java} -Xmx10g -jar \
//    ${params.gatk} \
//    -T IndelRealigner \
//    -R ${params.reference} \
//    -I split.bam \
//    -targetIntervals forIndelRealigner.intervals \
//    -o tmp2.bam

//  ${params.java} -Xmx10g -jar \
//    ${params.gatk} \
//    -T BaseRecalibrator \
//    -R ${params.reference} \
//    -I tmp2.bam \
//    -o recal_data.table \
//    -knownSites ${params.dbsnp} \
//    -nct 4

//  ${params.java} -Xmx10g -jar \
//    ${params.gatk} \
//    -T PrintReads \
//    -R ${params.reference} \
//    -I tmp2.bam \
//    --BQSR recal_data.table \
//    -o ${prefix}.gatk.bam
//
//  ${params.samtools} index ${prefix}.gatk.bam
//  """
//}

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
  clusterOptions '-l h_vmem=6G -pe smp 4'

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

//process cufflinks {
//  scratch true
//  executor 'sge'
//  clusterOptions '-l h_vmem=6G -pe smp 4'

// storeDir "${baseDir}/Cufflinks/${prefix}"

//  input:
//    set prefix, file(bam_file), file(bam_index_file) from gatk_out
//  output:
//    set prefix, file('*genes.fpkm_tracking'), file('*isoforms.fpkm_tracking'), file('*skipped.gtf'), file('*transcripts.gtf') into cufflinks_out
//    set file('*genes.fpkm_tracking') into cufflinks_out_genes
//    set file('*isoforms.fpkm_tracking') into cufflinks_out_isoforms

//  """
//  ${params.cufflinks} \
//    -q -p 4 \
//    -o . \
//    -b ${params.reference} \
//    -G ${params.gtf} \
//    ${bam_file}
//  mv genes.fpkm_tracking ${prefix}.genes.fpkm_tracking
//  mv isoforms.fpkm_tracking ${prefix}.isoforms.fpkm_tracking
//  mv skipped.gtf ${prefix}.skipped.gtf
//  mv transcripts.gtf ${prefix}.transcripts.gtf
//  """
//}

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

process convert_expression_files {

  storeDir "${baseDir}/CufflinksSTAR/"

  input:
    set genes from combine_cufflinks_genes

  output:
    set file('genes-symbols-fpkm.csv') into genes_symbols
    set file('genes-humanSymbols-fpkm.csv') into genes_symbols_human

  """
  Rscript ${params.seqpy}/ensembl2symbol_mouse.R \
    -data ${genes} \
    -out genes-symbols-fpkm.csv

  python ${params.seqpy}/mouseSymbol2human.py \
    --infile genes-symbols-fpkm.csv \
    --out genes-humanSymbols-fpkm.csv
  """
}

process htseq_reads {

  storeDir "${baseDir}/HTSeqCount-STAR/${prefix}"
  maxForks 16

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


process run_deseq2{

  storeDir "${baseDir}/DESeq2/"

  input:
    file counts_file from combine_htseq_out.first()
    set name, p1, p2, group1, group2 from deseq2_comparisons

  output:
    set file(name) into deseq2_out

  """
  python ${params.seqpy}/deseq2.py \
    --counts ${counts_file} \
    --group1 ${group1.join(",")} \
    --group2 ${group2.join(",")} \
    --phenotypes ${p1},${p2} \
    --outdir ${name}

  Rscript ${params.seqpy}/ensembl2symbol_mouse.R \
    -data ${name}/${name}_results.csv \
    -out ${name}/${name}_results_geneSymbols.csv
  """
}

process run_dexseq_count {
  storeDir "${baseDir}/DEXSeq/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_index_file) from star_out4
  output:
    set file('*.dexseq.counts') into dexseq_count_out
    set prefix into samples_combine_dexseq_count

  """
  python ${params.dexseq} \
    -p no \
    -s no \
    -f bam \
    ${params.dexseqgff} \
    ${bam_file} \
    ${prefix}.dexseq.counts
  """
}

process combine_dexseq_counts {
  storeDir "${baseDir}/DEXSeq/"
  input:
    file results from dexseq_count_out.toList()
    val all_samples from samples_combine_dexseq_count.toList()
  output:
    set file('DEXseq.counts.csv') into combine_dexseq_count_out

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out DEXseq.counts.csv
  """
}

process count_ptprk_rspo3_junction_reads {
  storeDir "${baseDir}/count_ptprk_rspo3/"

  input:
    file bams from bam_junction_count.toList()
    file bam_bai from bam_bai_junction_count.toList()
    file bamChimerics from star_chim_sam_out.toList()
    val all_samples from samples_junction_count.toList()

  output:
    set file('ptprk-rspo3_fusion-read-counts.csv') into fusion_out

  """
  python ${params.seqpy}/count_reads_for_junction.py \
    --samples ${all_samples.join(" ")} \
    --bam ${bams} \
    --bamChimeric ${bamChimerics} \
    --gene1chrom 10 \
    --gene1start 28075177 \
    --gene1end 28075178 \
    --gene2chrom 10 \
    --gene2start 29506579 \
    --gene2end 29506580 \
    --STAR ~/chm2059/lib/STAR-STAR_2.4.1d/bin/Linux_x86_64_static/STAR \
    --STARref ~/chm2059/dow_2016_5/result/RNAseq/STAR/ptprk-rspo3-fusion/star/ \
    --samtools ~/chm2059/lib/samtools-1.2/samtools \
    --out ptprk-rspo3_fusion-read-counts.csv \
    --gene1name Ptprk \
    --gene2name Rspo3 \
    --junctionPosition 97
  """
}

def getFastqFiles( Path directory ) {

  def prefix = directory.name

  def rl = '84'

  R1_files = []
  directory.eachFile{
    R1_files << it
  }
  R1_files.sort()
  return [prefix, R1_files, rl]
}
