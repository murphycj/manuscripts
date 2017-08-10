
basedir = "$baseDir"
params.seqpy = '/Users/charlesmurphy/Desktop/Research/lib/seqpy/bin'
params.count_file = '/Users/charlesmurphy/Desktop/Research/0816_sam/results/RNAseq_Samuel4427_2017_02_08/HTSeqCount/HTSeq.geneSymbols.counts.csv'


deseq2_comparisons = Channel
 .from(
   [
     'primary_vs_metastasis',
     'primary',
     'metastasis',
     ['7','8','9','10','19','20','21','22','26','27','35'],
     ['1','2','3','4','5','6','23','24','25','28','29','30','31','32','33','34','36','37','38','39','40','41','42','43','44','45','46']
   ],
   [
     'CINhighPrimary_vs_CINhighMetastasis',
     'CINhighPrimary',
     'CINhighMetastasis',
     ['7','26','27','35'],
     ['2','3','4','5','6','23','24','25','29','30','31','32']
   ],
   [
      'CIN_MedHigh_vs_STINGshRNA',
      'CIN_MedHigh',
      'STINGshRNA',
      ['13','14','17','18'],
      ['11','12']
   ],
   [
      'CIN_Med_vs_STINGshRNA',
      'CIN_Med',
      'STINGshRNA',
      ['17','18'],
      ['11','12']
   ],
   [
      'Control_vs_STINGshRNA',
      'Control',
      'STINGshRNA',
      ['17'],
      ['11','12']
   ],
   [
      'BoneMet_vs_BrainMet',
      'BoneMet',
      'BrainMet',
      ['4','6','28','29','33','39','43','44'],
      ['1','2','3','30','34','36','37','40','42','45','46']
   ],
   [
      'CINLowPrimary_vs_CINMedHighPrimary',
      'CINLowPrimary',
      'CINMedHighPrimary',
      ['8','9'],
      ['7','10','19','20','21','22','27','35']
   ],
   [
      'CINLowMets_vs_CINHighMets',
      'CINLowMets',
      'CINHighMets',
      ['34','39','45'],
      ['2','3','4','5','6','23','24','25','29','30','31','32']
   ]
 )

process run_deseq2{

  storeDir "${baseDir}"

  input:
    set name, p1, p2, group1, group2 from deseq2_comparisons

  output:
    file("${name}") into deseq2_out
    set name, file("${name}/*csv") into deseq2_out_csv

  """
  Rscript ${params.seqpy}/deseq2.R \
    -counts ${params.count_file} \
    -minCount 0 \
    -minRowCount 10 \
    -G1 ${group1.join(",")} \
    -G2 ${group2.join(",")} \
    -phenotype ${p1},${p2} \
    -outDir ${name}

  """
}

process run_rnk {

  storeDir "${baseDir}/${name}"

  input:
    set name, deseq2_out_file from deseq2_out_csv

  output:
    file("*.log2FC.rnk") into rnk_out_log2FC
    file("*statistic.rnk") into rnk_out_statistic
    file("*logPvalSignFC.rnk") into rnk_out_logPvalSignFC

  """
  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking log2FC \
    --out ${baseDir}/${name}/${name}.log2FC.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking statistic \
    --out ${baseDir}/${name}/${name}.statistic.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking logPvalSignFC \
    --out ${baseDir}/${name}/${name}.logPvalSignFC.rnk
  """
}
