
basedir = "$baseDir"
params.seqpy = '/Users/charlesmurphy/Desktop/Research/lib/seqpy/bin'
params.count_file = '/Users/charlesmurphy/Desktop/Research/0816_sam/results/RNAseq/HTSeqCount/HTSeq.gene.counts.csv'
params.count_file_symbol = '/Users/charlesmurphy/Desktop/Research/0816_sam/results/RNAseq/HTSeqCount/HTSeq.geneSymbols.counts.csv'

 deseq2_comparisons = Channel
  .from(
    [
      'CIN_low_vs_CIN_medium',
      'CIN_low',
      'CIN_medium',
      ['Kb','MK'],
      ['cont','Ka']
    ],
    [
      'CIN_low_vs_CIN_high',
      'CIN_low',
      'CIN_high',
      ['Kb','MK'],
      ['MKH']
    ],
    [
      'CIN_low_vs_CIN_rest',
      'CIN_low',
      'CIN_rest',
      ['Kb','MK'],
      ['MKH','cont','Ka']
    ],
  )


  deseq2_comparisons_symbols = Channel
   .from(
     [
       'CIN_low_vs_CIN_medium_symbols',
       'CIN_low',
       'CIN_medium',
       ['Kb','MK'],
       ['cont','Ka']
     ],
     [
       'CIN_low_vs_CIN_high_symbols',
       'CIN_low',
       'CIN_high',
       ['Kb','MK'],
       ['MKH']
     ],
     [
       'CIN_low_vs_CIN_rest_symbols',
       'CIN_low',
       'CIN_rest',
       ['Kb','MK'],
       ['MKH','cont','Ka']
     ],
   )


process run_deseq2{

  storeDir "${baseDir}/${name}"

  input:
    set name, p1, p2, group1, group2 from deseq2_comparisons

  output:
    set file('*pdf') into deseq2_out_pdf
    set name, file('*csv') into deseq2_out

  """
  python ${params.seqpy}/deseq2.py \
    --counts ${params.count_file} \
    --mincount 2 \
    --group1 ${group1.join(",")} \
    --group2 ${group2.join(",")} \
    --phenotypes ${p1},${p2} \
    --outdir ${name}

  """
}

process run_rnk {

  storeDir "${baseDir}/${name}"

  input:
    set name, file(deseq2_out_file) from deseq2_out

  output:
    set file("*.log2FC.rnk") into rnk_out_log2FC
    set file("*statistic.rnk") into rnk_out_statistic
    set file("*logPvalSignFC.rnk") into rnk_out_logPvalSignFC

  """
  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking log2FC \
    --out ${name}.log2FC.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking statistic \
    --out ${name}.statistic.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking logPvalSignFC \
    --out ${name}.logPvalSignFC.rnk
  """
}


process run_deseq2_symbols{

  storeDir "${baseDir}/${name}"

  input:
    set name, p1, p2, group1, group2 from deseq2_comparisons_symbols

  output:
    set file('*pdf') into deseq2_out_pdf2
    set name, file('*csv') into deseq2_out2

  """
  python ${params.seqpy}/deseq2.py \
    --counts ${params.count_file_symbol} \
    --mincount 2 \
    --group1 ${group1.join(",")} \
    --group2 ${group2.join(",")} \
    --phenotypes ${p1},${p2} \
    --outdir ${name}

  """
}

process run_rnk {

  storeDir "${baseDir}/${name}"

  input:
    set name, file(deseq2_out_file) from deseq2_out2

  output:
    set file("*.log2FC.rnk") into rnk_out_log2FC2
    set file("*statistic.rnk") into rnk_out_statistic2
    set file("*.abslog2FC.rnk") into rnk_out_log2FC2_abs
    set file("*.absstatistic.rnk") into rnk_out_statistic_abs
    set file("*logPvalSignFC.rnk") into rnk_out_logPvalSignFC2

  """
  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking log2FC \
    --out ${name}.log2FC.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking statistic \
    --out ${name}.statistic.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking log2FC \
    --abs \
    --out ${name}.abslog2FC.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking statistic \
    --abs \
    --out ${name}.absstatistic.rnk

  python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking logPvalSignFC \
    --out ${name}.logPvalSignFC.rnk
  """
}
