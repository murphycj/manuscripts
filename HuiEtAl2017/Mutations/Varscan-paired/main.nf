#!/usr/bin/env nextflow

basedir = "$baseDir"

params.bgzip = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/tabix-0.2.6/bgzip'
params.tabix = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/tabix-0.2.6/tabix'
params.bcftools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/bcftools-1.3.1/bcftools'
params.snpeff = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/snpEff'
params.vcftools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/vcftools_0.1.12b/bin/vcftools'
params.vcftools_bin = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/vcftools_0.1.12b/bin/'
params.java = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/jre1.8.0_25/bin/java'
params.varscan = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/VarScan.v2.3.7.jar'
params.gatk = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar'
params.seqpy = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/seqpy/bin/'
params.samtools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/samtools-1.3.1/samtools'

params.reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'
params.reference_dict = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.dict'
params.reference_fai = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa.fai'
params.dbsnp = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/mm10/mgp.v4.snps.dbSNP.vcf'
params.dbsnp_indels = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/mm10/mgp.v4.indels.dbSNP.vcf'
params.bed_file = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/data/Mouse_Exome_Design/110624_S0276129.bed'
params.samples = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/data/samples.xlsx'
params.rnaediting = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/RNAediting/gb-2012-13-4-r26-S1.mm10.sorted.vcf.gz'

bam_files = file("paired.csv")
pairs_varscan = Channel
  .from(bam_files)
  .splitCsv(header: false)
  .map{it -> [it[0], it[1], it[2], it[3], it[4], it[5]]}

process varscan {

  storeDir "${baseDir}/${tumor}"
  scratch true
  executor 'sge'
  stageInMode 'copy'
  clusterOptions '-l h_vmem=12G -pe smp 2 -l h_rt=48:00:00 -l athena=true'

  input:
    set tumor, control, tumor_bam, tumor_bam_bai, control_bam, control_bam_bai from pairs_varscan

  output:
    set tumor, control, file('*snps.vcf'), file('*indels.vcf') into varscan_out

  """
  rsync -L ${params.reference} ref.fa
  rsync -L ${params.reference_dict} ref.dict
  rsync -L ${params.reference_fai} ref.fa.fai
  rsync ${control_bam} control.bam
  rsync ${control_bam_bai} control.bam.bai
  rsync ${tumor_bam} tumor.bam
  rsync ${tumor_bam_bai} tumor.bam.bai

  ${params.samtools} mpileup \
    -d10000 \
    -Q 15 \
    -q 15 \
    -f ref.fa \
    -l ${params.bed_file} \
    control.bam \
    tumor.bam \
    | ${params.java} -jar ${params.varscan} somatic \
    --output-snp ${tumor}.snps.vcf \
    --output-indel ${tumor}.indels.vcf \
    --p-value 0.05 \
    --min-var-freq 0.05 \
    --min-coverage-tumor 12 \
    --min-coverage 12 \
    --output-vcf \
    --mpileup 1
  """
}

process varscan_process {
  storeDir "${baseDir}/${tumor}"

  maxForks 40

  input:
    set tumor, normal, vcf_snp, vcf_indel from varscan_out

  output:
    file('*somatic.hc.sorted.vcf.gz') into varscan_Somatic_hc_out
    file('*somatic.hc.sorted.vcf.gz.tbi') into varscan_Somatic_hc_out_index

    file('*somatic.sorted.vcf.gz') into varscan_Somatic_out
    file('*somatic.sorted.vcf.gz.tbi') into varscan_Somatic_out_index

    file('*.Germline.hc.sorted.vcf.gz') into varscan_process_germline
    file('*.Germline.hc.sorted.vcf.gz.tbi') into varscan_process_germline_index

    file('*.Germline.sorted.vcf.gz') into varscan_process_germline_hc
    file('*.Germline.sorted.vcf.gz.tbi') into varscan_process_germline_hc_index

  """
  rsync -L ${vcf_snp} ${tumor}.snps.vcf
  rsync -L ${vcf_indel} ${tumor}.indels.vcf

  ${params.vcftools} \
    --vcf ${tumor}.snps.vcf \
    --bed ${params.bed_file} \
    --recode --recode-INFO-all
  mv out.recode.vcf ${tumor}.snps.vcf
  ${params.vcftools} \
    --vcf ${tumor}.snps.vcf \
    --exclude-bed ${baseDir}/excluded_genes.sorted.bed \
    --recode --recode-INFO-all
  mv out.recode.vcf ${tumor}.snps.vcf

  ${params.vcftools} \
    --vcf ${tumor}.indels.vcf \
    --bed ${params.bed_file} \
    --recode --recode-INFO-all
  mv out.recode.vcf ${tumor}.indels.vcf
  ${params.vcftools} \
    --vcf ${tumor}.indels.vcf \
    --exclude-bed ${baseDir}/excluded_genes.sorted.bed \
    --recode --recode-INFO-all
  mv out.recode.vcf ${tumor}.indels.vcf

  echo ${normal} > name.txt
  echo ${tumor} >> name.txt

  ${params.bcftools} reheader ${tumor}.snps.vcf -s name.txt > tmp.vcf
  mv tmp.vcf ${tumor}.snps.vcf

  ${params.bcftools} reheader ${tumor}.indels.vcf -s name.txt > tmp.vcf
  mv tmp.vcf ${tumor}.indels.vcf

  java -jar ${params.varscan} somaticFilter \
    ${tumor}.snps.vcf \
    --indel-file ${tumor}.indels.vcf \
    --output-file ${tumor}.snps.filter.vcf \
    --p-value 0.01 \
    --min-coverage 12 \
    --min-reads2 4 \
    --min-var-freq 0.05
  java -jar ${params.varscan} somaticFilter ${tumor}.indels.vcf \
    --output-file ${tumor}.indels.filter.vcf \
    --p-value 0.01 \
    --min-coverage 12 \
    --min-reads2 4 \
    --min-var-freq 0.05

  java -jar ${params.varscan} processSomatic ${tumor}.snps.filter.vcf
  java -jar ${params.varscan} processSomatic ${tumor}.indels.filter.vcf

  ${params.vcftools} --vcf ${tumor}.snps.filter.Somatic.hc.vcf --remove-indv ${normal} --recode
  mv out.recode.vcf ${tumor}.snps.filter.Somatic.hc.tumor.vcf
  ${params.vcftools} --vcf ${tumor}.snps.filter.Somatic.vcf --remove-indv ${normal} --recode
  mv out.recode.vcf ${tumor}.snps.filter.Somatic.tumor.vcf
  ${params.vcftools} --vcf ${tumor}.snps.filter.Germline.hc.vcf --remove-indv ${tumor} --recode
  mv out.recode.vcf ${tumor}.snps.filter.Germline.hc.tumor.vcf
  ${params.vcftools} --vcf ${tumor}.snps.filter.Germline.vcf --remove-indv ${tumor} --recode
  mv out.recode.vcf ${tumor}.snps.filter.Germline.tumor.vcf

  ${params.vcftools} --vcf ${tumor}.indels.filter.Somatic.hc.vcf --remove-indv ${normal} --recode
  mv out.recode.vcf ${tumor}.indels.filter.Somatic.hc.tumor.vcf
  ${params.vcftools} --vcf ${tumor}.indels.filter.Somatic.vcf --remove-indv ${normal} --recode
  mv out.recode.vcf ${tumor}.indels.filter.Somatic.tumor.vcf
  ${params.vcftools} --vcf ${tumor}.indels.filter.Germline.hc.vcf --remove-indv ${tumor} --recode
  mv out.recode.vcf ${tumor}.indels.filter.Germline.hc.tumor.vcf
  ${params.vcftools} --vcf ${tumor}.indels.filter.Germline.vcf --remove-indv ${tumor} --recode
  mv out.recode.vcf ${tumor}.indels.filter.Germline.tumor.vcf

  ${params.vcftools_bin}/vcf-concat \
    ${tumor}.snps.filter.Somatic.hc.tumor.vcf \
    ${tumor}.indels.filter.Somatic.hc.tumor.vcf > ${tumor}.somatic.hc.vcf
  ${params.vcftools_bin}/vcf-sort ${tumor}.somatic.hc.vcf > ${tumor}.somatic.hc.sorted.vcf

  ${params.vcftools_bin}/vcf-concat \
    ${tumor}.snps.filter.Somatic.tumor.vcf \
    ${tumor}.indels.filter.Somatic.tumor.vcf > ${tumor}.somatic.vcf
  ${params.vcftools_bin}/vcf-sort ${tumor}.somatic.vcf > ${tumor}.somatic.sorted.vcf

  ${params.vcftools_bin}/vcf-concat \
    ${tumor}.snps.filter.Germline.hc.vcf \
    ${tumor}.indels.filter.Germline.hc.vcf > ${tumor}.Germline.hc.vcf
  ${params.vcftools_bin}/vcf-sort ${tumor}.Germline.hc.vcf > ${tumor}.Germline.hc.sorted.vcf

  ${params.vcftools_bin}/vcf-concat \
    ${tumor}.snps.filter.Germline.vcf \
    ${tumor}.indels.filter.Germline.vcf > ${tumor}.Germline.vcf
  ${params.vcftools_bin}/vcf-sort ${tumor}.Germline.vcf > ${tumor}.Germline.sorted.vcf

  ${params.bgzip} ${tumor}.somatic.hc.sorted.vcf
  ${params.tabix} ${tumor}.somatic.hc.sorted.vcf.gz

  ${params.bgzip} ${tumor}.somatic.sorted.vcf
  ${params.tabix} ${tumor}.somatic.sorted.vcf.gz

  ${params.bgzip} ${tumor}.Germline.hc.sorted.vcf
  ${params.tabix} ${tumor}.Germline.hc.sorted.vcf.gz

  ${params.bgzip} ${tumor}.Germline.sorted.vcf
  ${params.tabix} ${tumor}.Germline.sorted.vcf.gz

  """
}

somatic_hc = varscan_Somatic_hc_out.toList()
somatic_hc_index = varscan_Somatic_hc_out_index.toList()

germline = varscan_process_germline.toList()
germline_index = varscan_process_germline_index.toList()

germline_hc = varscan_process_germline_hc.toList()
germline_hc_index = varscan_process_germline_hc_index.toList()

process combine_vcf {

  storeDir "${baseDir}/"

  input:
    file somatic_hc
    file somatic_hc_index
    file germline
    file germline_index
    file germline_hc
    file germline_hc_index

  output:
    file('somatic_hc.varscan.vcf') into combine_vcf_out_somatic_hc
    set file('germline_hc.vcf.gz'), file('germline.vcf.gz') into combine_vcf_germline

  """
  export PATH=/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/tabix-0.2.6/:$PATH

  ${params.bcftools} merge -m none ${somatic_hc} > tmp.somatic.vcf
  ${params.bgzip} tmp.somatic.vcf
  ${params.tabix} tmp.somatic.vcf.gz
  ${params.vcftools_bin}/vcf-isec -c -f tmp.somatic.vcf.gz ${params.rnaediting} > somatic_hc.varscan.vcf
  rm tmp*
  ${params.bgzip} somatic_hc.varscan.vcf
  ${params.tabix} somatic_hc.varscan.vcf.gz

  ${params.bcftools} merge --force-samples -m none ${germline} > germline.vcf
  ${params.bgzip} germline.vcf
  ${params.tabix} germline.vcf.gz

  ${params.bcftools} merge --force-samples -m none ${germline_hc} > germline_hc.vcf
  ${params.bgzip} germline_hc.vcf
  ${params.tabix} germline_hc.vcf.gz

  ${params.vcftools_bin}/vcf-isec -c -f somatic_hc.varscan.vcf.gz germline.vcf.gz > somatic_hc.varscan.vcf
  """
}

process filter_somatic_snpeff {

  storeDir "${baseDir}/"

  input:
    file combine_vcf_out_somatic_hc

  output:
    file('somatic_hc.varscan.nodbsnp.vcf') into filter_somatic_snpeff_out
    file('somatic_hc.varscan.nodbsnp.vcf') into filter_somatic_snpeff_out2
    file('somatic_hc.varscan.nodbsnp.vcf') into filter_somatic_snpeff_out3

  """
  ${params.java} -Xmx40g -jar \
    ${params.snpeff}/SnpSift.jar annotate \
    ${params.dbsnp} \
    ${combine_vcf_out_somatic_hc} \
    | ${params.java} -Xmx40g -jar \
    ${params.snpeff}/SnpSift.jar annotate \
    ${params.dbsnp_indels} \
    | ${params.java} -Xmx40g -jar \
    ${params.snpeff}/snpEff.jar eff \
    -no-downstream \
    -no-upstream \
    -v GRCm38.86 \
    -canon \
    | ${params.java} -Xmx40g -jar \
    ${params.snpeff}/SnpSift.jar filter \
    "! exists ID" > somatic_hc.varscan.nodbsnp.vcf
  """
}

process vcf_to_bed {

  storeDir "${baseDir}/"

  input:
    file filter_somatic_snpeff_out

  output:
    file('variant_positions.bed') into variant_positions
    file('variant_positions.bed') into variant_positions2

  """
  export PATH=/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/bedops_linux_x86_64-v2.4.2/:$PATH
  /athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/bedops_linux_x86_64-v2.4.2/vcf2bed \
    --do-not-sort < ${filter_somatic_snpeff_out} \
    | awk -v OFS='\t' '{print \$1, \$2, \$3}' > variant_positions.bed
  """
}

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

process get_pileup {

  storeDir "${baseDir}/mpileups/${s}"
  scratch true
  executor 'sge'
  stageInMode 'copy'
  clusterOptions '-l h_vmem=4G -pe smp 2 -l h_rt=1:00:00 -l athena=true'

  input:
    set s, bam, bam_bai from bams
    file variant_positions

  output:
    val s into all_pileup_samples
    file('*.pileup') into all_pileup_files

  """
  rsync -L ${params.reference} ref.fa
  rsync -L ${params.reference_dict} ref.dict
  rsync -L ${params.reference_fai} ref.fa.fai
  ${params.samtools} mpileup \
    -f ref.fa \
    -l ${variant_positions} \
    -Q 10 \
    -q 10 \
    ${bam} \
    -o ${s}.pileup
  """
}

process combine_pileup {

  storeDir "${baseDir}"

  input:
    file ap from all_pileup_files.toList()
    val s from all_pileup_samples.toList()
    file variant_positions2
    file filter_somatic_snpeff_out2

  output:
    file('coverage.csv') into combine_pileup_coverage
    file('support.csv') into combine_pileup_support

  """
  python ${params.seqpy}/aggregate_mpileup.py \
    --vcf ${filter_somatic_snpeff_out2} \
    --pileups ${ap} \
    --samples ${s.join(' ')}
  """
}


process filter {

  storeDir "${baseDir}"

  input:
    file filter_somatic_snpeff_out3
    file combine_pileup_coverage
    file combine_pileup_support

  output:
    file('somatic.nodbsnp.varscan.filtered.vcf') into filter_out
    file('somatic.nodbsnp.varscan.removed.vcf') into filter_out_removed

  """
  python ${baseDir}/filter_vcf.py \
    --vcf ${filter_somatic_snpeff_out3} \
    --samples ${params.samples} \
    --Tn 212 \
    --Tr 4 \
    --Cn 10 \
    --Cc 12 \
    --Cr 4 \
    --coverage ${combine_pileup_coverage} \
    --support ${combine_pileup_support} \
    --out somatic.nodbsnp.varscan.filtered.vcf
  """
}

process vcf_to_table {

  storeDir "${baseDir}"

  input:
    file filter_out

  output:
    file('somatic.nodbsnp.varscan.filtered.csv') into vcf_to_table_out
    file('somatic.nodbsnp.varscan.filtered.allEffects.csv') into vcf_to_table_out_all

  """
  python ${params.seqpy}/vcf_to_table.py \
    --vcf ${filter_out} \
    --out somatic.nodbsnp.varscan.filtered.csv \
    --varscan

  python ${params.seqpy}/vcf_to_table.py \
    --vcf ${filter_out} \
    --everything \
    --out somatic.nodbsnp.varscan.filtered.allEffects.csv \
    --varscan
  """
}
