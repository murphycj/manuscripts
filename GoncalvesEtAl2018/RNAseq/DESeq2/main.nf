
basedir = "$baseDir"
params.seqpy = '/Users/charlesmurphy/Desktop/Research/lib/seqpy/bin'
params.count_file = '/Users/charlesmurphy/Desktop/Research/0815_jihye/results/RNAseq/HTSeqCount/rename_filtered/HTSeq.human.geneSymbols.counts.csv'

deseq2_comparisons = Channel
  .from(
    [
      'Liver_C_vs_F',
      'C',
      'F',
      ['C_L1','C_L2','C_L3','C_L4'],
      ['F_L1','F_L2','F_L3','F_L4']
    ],
    [
      'Liver_C_vs_G',
      'C',
      'G',
      ['C_L1','C_L2','C_L3','C_L4'],
      ['G_L1','G_L2','G_L3','G_L4']
    ],
    [
      'Liver_C_vs_H',
      'C',
      'H',
      ['C_L1','C_L2','C_L3','C_L4'],
      ['H_L1','H_L2','H_L3','H_L4']
    ],
    [
      'Liver_C_vs_W',
      'C',
      'W',
      ['C_L1','C_L2','C_L3','C_L4'],
      ['W_L1','W_L2','W_L3','W_L4']
    ],
    [
      'Liver_F_vs_G',
      'F',
      'G',
      ['F_L1','F_L2','F_L3','F_L4'],
      ['G_L1','G_L2','G_L3','G_L4']
    ],
    [
      'Liver_F_vs_H',
      'F',
      'H',
      ['F_L1','F_L2','F_L3','F_L4'],
      ['H_L1','H_L2','H_L3','H_L4']
    ],
    [
      'Liver_G_vs_H',
      'G',
      'H',
      ['G_L1','G_L2','G_L3','G_L4'],
      ['H_L1','H_L2','H_L3','H_L4']
    ],
    [
      'Liver_W_vs_H',
      'W',
      'H',
      ['W_L1','W_L2','W_L3','W_L4'],
      ['H_L1','H_L2','H_L3','H_L4']
    ],
    [
      'Liver_GF_vs_H',
      'GF',
      'H',
      ['F_L1','F_L2','F_L3','F_L4','G_L1','G_L2','G_L3','G_L4'],
      ['H_L1','H_L2','H_L3','H_L4']
    ],
    [
      'Liver_CGF_vs_H',
      'CGF',
      'H',
      ['C_L1','C_L2','C_L3','C_L4','F_L1','F_L2','F_L3','F_L4','G_L1','G_L2','G_L3','G_L4'],
      ['H_L1','H_L2','H_L3','H_L4']
    ],
    [
      'Liver_CGFW_vs_H',
      'CGFW',
      'H',
      ['C_L1','C_L2','C_L3','C_L4','F_L1','F_L2','F_L3','F_L4','G_L1','G_L2','G_L3','G_L4','W_L1','W_L2','W_L3','W_L4'],
      ['H_L1','H_L2','H_L3','H_L4']
    ],
    [
      'Liver_GF_vs_W',
      'GF',
      'W',
      ['F_L1','F_L2','F_L3','F_L4','G_L1','G_L2','G_L3','G_L4'],
      ['W_L1','W_L2','W_L3','W_L4']
    ],
    [
      'Liver_CGF_vs_W',
      'CGF',
      'W',
      ['F_L1','F_L2','F_L3','F_L4','G_L1','G_L2','G_L3','G_L4','C_L1','C_L2','C_L3','C_L4'],
      ['W_L1','W_L2','W_L3','W_L4']
    ],
    [
      'Liver_CGFH_vs_W',
      'CGFH',
      'W',
      ['F_L1','F_L2','F_L3','F_L4','G_L1','G_L2','G_L3','G_L4','C_L1','C_L2','C_L3','C_L4','H_L1','H_L2','H_L3','H_L4'],
      ['W_L1','W_L2','W_L3','W_L4']
    ],
    [
      'Liver_GH_vs_F',
      'GH',
      'F',
      ['G_L1','G_L2','G_L3','G_L4','H_L1','H_L2','H_L3','H_L4'],
      ['F_L1','F_L2','F_L3','F_L4']
    ],
    [
      'Liver_CGH_vs_F',
      'CGH',
      'F',
      ['G_L1','G_L2','G_L3','G_L4','H_L1','H_L2','H_L3','H_L4','C_L1','C_L2','C_L3','C_L4'],
      ['F_L1','F_L2','F_L3','F_L4']
    ],
    [
      'Liver_FH_vs_G',
      'FH',
      'G',
      ['F_L1','F_L2','F_L3','F_L4','H_L1','H_L2','H_L3','H_L4'],
      ['G_L1','G_L2','G_L3','G_L4']
    ],
    [
      'Liver_CFH_vs_G',
      'CFH',
      'G',
      ['F_L1','F_L2','F_L3','F_L4','H_L1','H_L2','H_L3','H_L4','C_L1','C_L2','C_L3','C_L4'],
      ['G_L1','G_L2','G_L3','G_L4']
    ],
    [
      'Tumor_C_vs_F',
      'C',
      'F',
      ['C1','C2','C3','C4'],
      ['F1','F2','F3','F4']
    ],
    [
      'Tumor_C_vs_G',
      'C',
      'G',
      ['C1','C2','C3','C4'],
      ['G1','G2','G3','G4']
    ],
    [
      'Tumor_C_vs_H',
      'C',
      'H',
      ['C1','C2','C3','C4'],
      ['H1','H2','H3','H4']
    ],
    [
      'Tumor_C_vs_W',
      'C',
      'W',
      ['C1','C2','C3','C4'],
      ['W1','W2','W3','W4']
    ],
    [
      'Tumor_C_vs_W_noW3',
      'C',
      'W',
      ['C1','C2','C3','C4'],
      ['W1','W2','W4']
    ],
    [
      'Tumor_F_vs_G',
      'F',
      'G',
      ['F1','F2','F3','F4'],
      ['G1','G2','G3','G4']
    ],
    [
      'Tumor_F_vs_H',
      'F',
      'H',
      ['F1','F2','F3','F4'],
      ['H1','H2','H3','H4']
    ],
    [
      'Tumor_G_vs_H',
      'G',
      'H',
      ['G1','G2','G3','G4'],
      ['H1','H2','H3','H4']
    ],
    [
      'Tumor_GF_vs_H',
      'G_and_F',
      'H',
      ['G1','G2','G3','G4','F1','F2','F3','F4'],
      ['H1','H2','H3','H4']
    ],
    [
      'Tumor_H_vs_W',
      'H',
      'W',
      ['H1','H2','H3','H4'],
      ['W1','W2','W4']
    ],
    [
      'Tumor_FG_vs_W',
      'F_and_G',
      'W',
      ['G1','G2','G3','G4','F1','F2','F3','F4'],
      ['W1','W2','W4']
    ],
    [
      'Tumor_HFG_vs_W',
      'HFG',
      'W',
      ['G1','G2','G3','G4','F1','F2','F3','F4','H1','H2','H3','H4'],
      ['W1','W2','W4']
    ],
    [
      'Tumor_WFG_vs_H',
      'WFG',
      'H',
      ['G1','G2','G3','G4','F1','F2','F3','F4','W1','W2','W4'],
      ['H1','H2','H3','H4']
    ],
    [
      'Tumor_CFG_vs_H',
      'CFG',
      'H',
      ['C1','C2','C3','C4','G1','G2','G3','G4','F1','F2','F3','F4'],
      ['H1','H2','H3','H4']
    ],
    [
      'Tumor_CFG_vs_W',
      'CFG',
      'W',
      ['C1','C2','C3','C4','G1','G2','G3','G4','F1','F2','F3','F4'],
      ['W1','W2','W4']
    ],
    [
      'Tumor_CFGH_vs_W',
      'CFGH',
      'W',
      ['C1','C2','C3','C4','G1','G2','G3','G4','F1','F2','F3','F4','H1','H2','H3','H4'],
      ['W1','W2','W4']
    ],
    [
      'Tumor_CFGW_vs_H',
      'CFGW',
      'H',
      ['C1','C2','C3','C4','G1','G2','G3','G4','F1','F2','F3','F4','W1','W2','W4'],
      ['H1','H2','H3','H4']
    ],
    [
      'Colon_C_vs_F',
      'C',
      'F',
      ['C1-Colon','C2-Colon','C3-Colon','C4-Colon'],
      ['F1-Colon','F2-Colon','F3-Colon','F4-Colon']
    ],
    [
      'Colon_C_vs_G',
      'C',
      'G',
      ['C1-Colon','C2-Colon','C3-Colon','C4-Colon'],
      ['G1-Colon','G2-Colon','G3-Colon']
    ],
    [
      'Colon_C_vs_H',
      'C',
      'H',
      ['C1-Colon','C2-Colon','C3-Colon','C4-Colon'],
      ['H1-Colon','H2-Colon','H3-Colon','H4-Colon']
    ],
    [
      'Colon_F_vs_G',
      'F',
      'G',
      ['F1-Colon','F2-Colon','F3-Colon','F4-Colon'],
      ['G1-Colon','G2-Colon','G3-Colon']
    ],
    [
      'Colon_F_vs_H',
      'F',
      'H',
      ['F1-Colon','F2-Colon','F3-Colon','F4-Colon'],
      ['H1-Colon','H2-Colon','H3-Colon','H4-Colon']
    ],
    [
      'Colon_G_vs_H',
      'G',
      'H',
      ['G1-Colon','G2-Colon','G3-Colon'],
      ['H1-Colon','H2-Colon','H3-Colon','H4-Colon']
    ],
    [
      'Colon_GF_vs_H',
      'GF',
      'H',
      ['G1-Colon','G2-Colon','G3-Colon','F1-Colon','F2-Colon','F3-Colon','F4-Colon'],
      ['H1-Colon','H2-Colon','H3-Colon','H4-Colon']
    ],
    [
      'Colon_GFC_vs_H',
      'GFC',
      'H',
      ['G1-Colon','G2-Colon','G3-Colon','F1-Colon','F2-Colon','F3-Colon','F4-Colon','C1-Colon','C2-Colon','C3-Colon','C4-Colon'],
      ['H1-Colon','H2-Colon','H3-Colon','H4-Colon']
    ],
    [
      'Colon_CGH_vs_F',
      'CGH',
      'F',
      ['G1-Colon','G2-Colon','G3-Colon','C1-Colon','C2-Colon','C3-Colon','C4-Colon','H1-Colon','H2-Colon','H3-Colon','H4-Colon'],
      ['F1-Colon','F2-Colon','F3-Colon','F4-Colon']
    ],
    [
      'Colon_FH_vs_G',
      'FH',
      'G',
      ['F1-Colon','F2-Colon','F3-Colon','F4-Colon','H1-Colon','H2-Colon','H3-Colon','H4-Colon'],
      ['G1-Colon','G2-Colon','G3-Colon']
    ],
    [
      'Colon_CFH_vs_G',
      'CFH',
      'G',
      ['F1-Colon','F2-Colon','F3-Colon','F4-Colon','H1-Colon','H2-Colon','H3-Colon','H4-Colon','C1-Colon','C2-Colon','C3-Colon','C4-Colon'],
      ['G1-Colon','G2-Colon','G3-Colon']
    ]
  )

process run_deseq2{

  storeDir "${baseDir}"

  input:
    set name, p1, p2, group1, group2 from deseq2_comparisons

  output:
    file("${name}") into deseq2_out_pdf
    set name, file("${name}/*csv") into deseq2_out

  """
  Rscript ${params.seqpy}/deseq2.R \
    -counts ${params.count_file} \
    -G1 ${group1.join(",")} \
    -G2 ${group2.join(",")} \
    -minCount 2 \
    -minRowCount 2 \
    -phenotype ${p1},${p2} \
    -outDir ${name}

  """
}


process run_rnk {

  storeDir "${baseDir}/${name}"

  input:
    set name, file(deseq2_out_file) from deseq2_out

  output:
    file("*.log2FC.rnk") into rnk_out_log2FC
    file("*statistic.rnk") into rnk_out_statistic
    file("*logPvalSignFC.rnk") into rnk_out_logPvalSignFC

  """
  python ${params.seqpy}/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking log2FC \
    --out ${name}.log2FC.rnk

  python ${params.seqpy}/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking statistic \
    --out ${name}.statistic.rnk

  python ${params.seqpy}/deseq2rnk.py \
    --infile ${deseq2_out_file} \
    --ranking logPvalSignFC \
    --out ${name}.logPvalSignFC.rnk
  """
}
