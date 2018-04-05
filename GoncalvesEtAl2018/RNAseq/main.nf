#!/usr/bin/env nextflow

PE50 = (1..31).collect{
  'M' + it
}
PE75 = (1..20).collect{
  'C' + it
}
human = (1..20).collect{'A' + it}  + (1..20).collect{'B' + it}


params.fastqs = "/home/chm2059/chm2059/jihye_2015_8/data/RNAseq/*"

basedir = "$baseDir"

params.mouse_reference = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.mouse_reference_dict = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.dna.primary_assembly.dict'
params.mouse_reference_fai = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'

params.human_reference = '/home/chm2059/chm2059/data/refdata/GRCh38_91/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
params.human_reference_dict = '/home/chm2059/chm2059/data/refdata/GRCh38_91/Homo_sapiens.GRCh38.dna.primary_assembly.dict'
params.human_reference_fai = '/home/chm2059/chm2059/data/refdata/GRCh38_91/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai'

params.dexseq_gff_mouse = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.91.chr.dexseq.gtf'
params.dexseq_gff_human = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCh38_91/Homo_sapiens.GRCh38.91.chr.dexseq.gtf'

params.star = '/home/chm2059/chm2059/lib/STAR-STAR_2.4.1d/bin/Linux_x86_64_static/STAR'
params.outFilterMismatchNmax = 5
params.star_threads = 8

fqfiles_fastqc = Channel.create()
fqfiles_star = Channel.create()

Channel
  .fromPath(params.fastqs, type: 'dir')
  .map { path ->
       (prefix, R1_files, R2_files, star, pe, overhang, reference, gtf) = getFastqFiles(path)
       tuple(prefix, R1_files, R2_files, star, pe, overhang, reference, gtf)
  }
  .filter{
    it[0]!='Summary'
  }
  .filter{
    it[1].size()>0
  }
  .separate(fqfiles_fastqc, fqfiles_star){it->[it,it]}

process fastqc {

  storeDir "${baseDir}/fastqc/${prefix}"
  maxForks 6

  input:
    set prefix, file(read1), file(read2), star, pe, overhang, reference, gtf from fqfiles_fastqc

  output:
    file('*zip') into fastqc_results
    val prefix into samples

  script:
  if ( pe == 'no' )
    """
    zcat ${read1} | ${params.fastqc} stdin --outdir=.
    mv stdin_fastqc.zip ${prefix}.zip
    """
  else
    """
    zcat ${read1} ${read2} | ${params.fastqc} stdin --outdir=.
    mv stdin_fastqc.zip ${prefix}.zip
    """
}


process star {
  maxForks 4

  storeDir "${baseDir}/STAR/${prefix}"

  input:
    set prefix, file(read1), file(read2), star, pe, overhang, reference, gtf from fqfiles_star

  output:
    set prefix, file('*bam'), file('*bam.bai'), reference, gtf into star_out
    set prefix, file('*bam'), file('*bam.bai'), reference, gtf into star_out2
    set prefix, file('*bam'), file('*bam.bai'), reference, gtf into star_out3
    set prefix, file('*bam'), file('*bam.bai'), reference, gtf into star_out_varscan_snp
    set prefix, file('*bam'), file('*bam.bai'), reference, gtf into star_out_varscan_indel
    set prefix, file('*bam'), file('*bam.bai'), reference, gtf into star_out_gatk
    set prefix, file('*bam'), file('*bam.bai'), reference, gtf, pe into star_out_dexseq
    file('*Log.final.out') into star_log
    val prefix into samples_star_qc
    file('*SJ.out.tab') into star_out_SJ
    file('*Log.progress.out') into star_out_progress


  script:
  read1=read1.toString().replaceAll(/ /,",")
  read2=read2.toString().replaceAll(/ /,",")

  if (pe=='yes')
    """

    mkdir 1pass
    mkdir star_2pass
    ${params.star} \
      --outSAMtype BAM SortedByCoordinate \
      --runThreadN ${params.star_threads} \
      --outSAMstrandField intronMotif \
      --outFilterMismatchNmax ${params.outFilterMismatchNmax} \
      --genomeDir ${star} \
      --readFilesCommand zcat \
      --readFilesIn ${read1} ${read2} \
      --outFileNamePrefix ./1pass/tmp.

    ${params.star} \
      --runMode genomeGenerate \
      --genomeDir ./star_2pass \
      --genomeFastaFiles ${reference} \
      --sjdbFileChrStartEnd ./1pass/tmp.SJ.out.tab \
      --sjdbOverhang ${overhang} \
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
  else
    """

    mkdir 1pass
    mkdir star_2pass
    ${params.star} \
      --outSAMtype BAM SortedByCoordinate \
      --runThreadN ${params.star_threads} \
      --outSAMstrandField intronMotif \
      --outFilterMismatchNmax ${params.outFilterMismatchNmax} \
      --genomeDir ${star} \
      --readFilesCommand zcat \
      --readFilesIn ${read1} \
      --outFileNamePrefix ./1pass/tmp.

    ${params.star} \
      --runMode genomeGenerate \
      --genomeDir ./star_2pass \
      --genomeFastaFiles ${reference} \
      --sjdbFileChrStartEnd ./1pass/tmp.SJ.out.tab \
      --sjdbOverhang ${overhang} \
      --runThreadN ${params.star_threads}

    ${params.star} \
      --outSAMtype BAM SortedByCoordinate \
      --runThreadN ${params.star_threads} \
      --outSAMstrandField intronMotif \
      --genomeDir ./star_2pass \
      --readFilesCommand zcat \
      --readFilesIn ${read1} \
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
    set prefix, file(bam_file), file(bam_bai_file), reference, gtf from star_out3

  output:
    file('*flagstat') into flagstat_out

  """
  ${params.samtools} flagstat ${bam_file} > ${prefix}.flagstat

  """
}

process cufflinks_star {
  scratch true
  executor 'sge'
  clusterOptions '-l h_vmem=3G -pe smp 6 -l athena=true'

  storeDir "${baseDir}/CufflinksSTAR/${prefix}"

  input:
    set prefix, file(bam_file), file(bam_bai_file), reference, gtf from star_out2
  output:
    set prefix, file('*genes.fpkm_tracking'), file('*isoforms.fpkm_tracking'), file('*skipped.gtf'), file('*transcripts.gtf') into cufflinks_out

    """
    ${params.cufflinks} \
      -q -p 6 \
      -o . \
      -b ${reference} \
      -G ${gtf} \
      ${bam_file}
    mv genes.fpkm_tracking ${prefix}.genes.fpkm_tracking
    mv isoforms.fpkm_tracking ${prefix}.isoforms.fpkm_tracking
    mv skipped.gtf ${prefix}.skipped.gtf
    mv transcripts.gtf ${prefix}.transcripts.gtf
    """
}

cufflinks_out_human = Channel.create()
cufflinks_out_mouse = Channel.create()

cufflinks_out.choice(cufflinks_out_human,cufflinks_out_mouse) { a -> (a[0] =~ /^A/) || (a[0] =~ /^B/) ? 0 : 1 }

cufflinks_out_genes_human = Channel.create()
cufflinks_out_isoforms_human = Channel.create()
samples_cufflinks_human = Channel.create()

cufflinks_out_genes_mouse = Channel.create()
cufflinks_out_isoforms_mouse = Channel.create()
samples_cufflinks_mouse = Channel.create()

cufflinks_out_human
  .separate(samples_cufflinks_human,cufflinks_out_genes_human, cufflinks_out_isoforms_human) {a -> [a[0],a[1],a[2]]}

cufflinks_out_mouse
  .separate(samples_cufflinks_mouse,cufflinks_out_genes_mouse,cufflinks_out_isoforms_mouse) {a -> [a[0],a[1],a[2]]}

process combine_cufflinks_human {

  storeDir "${baseDir}/CufflinksSTAR/"

  input:
    file genes from cufflinks_out_genes_human.toList()
    file isoforms from cufflinks_out_isoforms_human.toList()
    val all_samples from samples_cufflinks_human.toList()

  output:
    set file('isoforms-humanFPKM.csv'), file('genes-humanFPKM.csv') into combine_cufflinks_genes_human

  """
  python /home/chm2059/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${genes} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out genes-humanFPKM.csv

  python /home/chm2059/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${isoforms} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out isoforms-humanFPKM.csv
  """
}

process combine_cufflinks_mouse {

  storeDir "${baseDir}/CufflinksSTAR/"

  input:
    file genes from cufflinks_out_genes_mouse.toList()
    file isoforms from cufflinks_out_isoforms_mouse.toList()
    val all_samples from samples_cufflinks_mouse.toList()

  output:
    set file('isoforms-mouseFPKM.csv'), file('genes-mouseFPKM.csv') into combine_cufflinks_genes_mouse

  """
  python /home/chm2059/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${genes} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out genes-mouseFPKM.csv

  python /home/chm2059/chm2059/lib/seqpy/bin/aggregate_fpkm.py \
    --files ${isoforms} \
    --samples ${all_samples.join(" ")} \
    --duplicates highest \
    --out isoforms-mouseFPKM.csv
  """
}


process htseq_reads {

  storeDir "${baseDir}/HTSeqCount/${prefix}"
  maxForks 20

  input:
    set prefix, file(bam_file), file(bam_index_file), reference, gtf from star_out

  output:
    set prefix, file('*count') into htseq_reads_out

  """
  python ${params.htseq} \
    -s no \
    -f bam \
    ${bam_file} \
    ${gtf} > ${prefix}.count
  """
}

htseq_out_human = Channel.create()
htseq_out_mouse = Channel.create()

htseq_reads_out.choice(htseq_out_human,htseq_out_mouse) { a -> (a[0] =~ /^A/) || (a[0] =~ /^B/) ? 0 : 1 }

htseq_out_genes_human = Channel.create()
samples_htseq_human = Channel.create()

htseq_out_genes_mouse = Channel.create()
samples_htseq_mouse = Channel.create()

htseq_out_human
  .separate(samples_htseq_human, htseq_out_genes_human) {a -> [a[0],a[1]]}

htseq_out_mouse
  .separate(samples_htseq_mouse, htseq_out_genes_mouse) {a -> [a[0],a[1]]}

process combine_htseq_human {

  storeDir "${baseDir}/HTSeqCount/"

  input:
    file results from htseq_out_genes_human.toList()
    val all_samples from samples_htseq_human.toList()

  output:
    file('HTSeq.gene.humanCounts.csv') into combine_htseq_out_human_count

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out HTSeq.gene.humanCounts.csv
  """
}

process combine_htseq_mouse {

  storeDir "${baseDir}/HTSeqCount/"

  input:
    file results from htseq_out_genes_mouse.toList()
    val all_samples from samples_htseq_mouse.toList()

  output:
    file('HTSeq.gene.mouseCounts.csv') into combine_htseq_out_mouse_count

  """
  python ${params.seqpy}/aggregate_htseqcount.py \
    --files ${results} \
    --samples ${all_samples.join(" ")} \
    --out HTSeq.gene.mouseCounts.csv
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

def getFastqFiles( Path directory ) {

  def prefix = directory.name

  def pe = ''
  def rl = ''

  if (prefix in PE50) {
    star = '/home/chm2059/chm2059/data/refdata/GRCm38_91/star100'
    pe = 'yes'
    overhang = 50
    reference = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.dna.primary_assembly.fa'
    gtf = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.91.chr.gtf'
  } else if (prefix in PE75){
    star = '/home/chm2059/chm2059/data/refdata/GRCm38_91/star100'
    pe = 'yes'
    overhang = 74
    reference = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.dna.primary_assembly.fa'
    gtf = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.91.chr.gtf'
  } else if (prefix in human) {
    star = '/home/chm2059/chm2059/data/refdata/GRCh38/star84'
    pe = 'no'
    overhang = 84
    reference = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    gtf = '/home/chm2059/chm2059/data/refdata/GRCh38/Homo_sapiens.GRCh38.79.gtf'
  } else {
    star = '/home/chm2059/chm2059/data/refdata/GRCm38_91/star100'
    pe = 'no'
    overhang = 84
    reference = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.dna.primary_assembly.fa'
    gtf = '/home/chm2059/chm2059/data/refdata/GRCm38_91/Mus_musculus.GRCm38.91.chr.gtf'
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
  return [prefix, R1_files, R2_files, star, pe, overhang, reference, gtf]
}
