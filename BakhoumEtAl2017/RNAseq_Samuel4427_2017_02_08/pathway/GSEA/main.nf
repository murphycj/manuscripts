
params.seqpy = '/home/chm2059/chm2059/lib/seqpy/bin/'
params.gsea = '/home/chm2059/chm2059/lib/gsea2-2.2.1.jar'

msigdb = Channel
  .from(
    'c2.cp.v6.0.symbols.gmt',
    'h.all.v6.0.symbols.gmt',
    'c1.all.v6.0.symbols.gmt'
  )

gsea_comparisons = Channel
 .from(
   [
    'CIN_MedHigh_vs_STINGshRNA',
    'log2FC'
   ],
   [
    'CINLowMets_vs_CINHighMets',
    'log2FC'
   ],
   [
    'CINLowPrimary_vs_CINMedHighPrimary',
    'log2FC'
   ],
   [
    'BoneMet_vs_BrainMet',
    'log2FC',
   ],
   [
    'CINhighPrimary_vs_CINhighMetastasis',
    'log2FC'
   ],
   [
    'Control_vs_STINGshRNA',
    'log2FC'
   ],
   [
    'primary_vs_metastasis',
    'log2FC'
   ]
 )
 .spread(msigdb)

process gsea {

  storeDir "${comparison}/${rnk}/"
  maxForks 15

  input:
    set comparison, rnk, gmt from gsea_comparisons

  output:
    set comparison, rnk, gmt, file("${comparison}_${rnk}_${gmt}*") into gsea_out

  """
  java -cp ${params.gsea} -Xmx2048m xtools.gsea.GseaPreranked \
    -rnk /home/chm2059/chm2059/0816_sam/results/RNAseq_Samuel4427_2017_02_08/DESeq2/${comparison}/${comparison}.${rnk}.rnk \
    -gmx /home/chm2059/chm2059/data/msigdb/msigdb_v6.0_GMTs/${gmt} \
    -collapse false \
    -mode Max_probe \
    -norm meandiv \
    -nperm 1000 \
    -scoring_scheme weighted \
    -rpt_label ${comparison}_${rnk}_${gmt} \
    -include_only_symbols true \
    -make_sets true \
    -plot_top_x 40 \
    -rnd_seed timestamp \
    -set_max 500 \
    -set_min 15 \
    -zip_report true \
    -out . \
    -gui false

  """
}

process get_output {

  storeDir "${baseDir}/${comparison}"
  maxForks 15

  input:
    set comparison, rnk, gmt, fin from gsea_out

  output:
    file("${comparison}_${rnk}_${gmt}.xlsx") into get_output_out

  """
  python ${params.seqpy}/get_gsea_output.py \
    --indir ${fin} \
    --up_tab_name up_regulated \
    --down_tab_name down_regulated \
    --allGenes \
    --leadingEdgeGenes \
    --leadingEdge \
    --out ${comparison}_${rnk}_${gmt}
  """
}
