
params.seqpy = '/home/chm2059/chm2059/lib/seqpy/bin/'
params.gsea = '/home/chm2059/chm2059/lib/gsea2-2.2.1.jar'

msigdb = Channel
  .from(
    'c2.cp.v5.1.symbols.gmt',
    'c3.tft.v5.1.symbols.gmt',
    'c5.all.v5.1.symbols.gmt',
    'c2.all.v5.1.symbols.gmt',
    'c6.all.v5.1.symbols.gmt',
    'h.all.v5.1.symbols.gmt',
    'c7.all.v5.1.symbols.gmt'
  )

gsea_comparisons = Channel
 .from(
   [
     'CIN_low_vs_CIN_medium_symbols',
     'log2FC'
   ],
   [
     'CIN_low_vs_CIN_medium_symbols',
     'statistic'
   ],
   [
     'CIN_low_vs_CIN_medium_symbols',
     'logPvalSignFC'
   ],
   [
     'CIN_low_vs_CIN_medium_symbols',
     'absstatistic'
   ],
   [
     'CIN_low_vs_CIN_medium_symbols',
     'abslog2FC'
   ],
   [
     'CIN_low_vs_CIN_high_symbols',
     'log2FC'
   ],
   [
     'CIN_low_vs_CIN_high_symbols',
     'statistic'
   ],
   [
     'CIN_low_vs_CIN_high_symbols',
     'abslog2FC'
   ],
   [
     'CIN_low_vs_CIN_high_symbols',
     'absstatistic'
   ],
   [
     'CIN_low_vs_CIN_high_symbols',
     'logPvalSignFC'
   ],
   [
     'CIN_low_vs_CIN_rest_symbols',
     'statistic'
   ],
   [
     'CIN_low_vs_CIN_rest_symbols',
     'log2FC'
   ],
   [
     'CIN_low_vs_CIN_rest_symbols',
     'absstatistic'
   ],
   [
     'CIN_low_vs_CIN_rest_symbols',
     'abslog2FC'
   ],
   [
    'CIN_low_vs_CIN_rest_symbols',
    'logPvalSignFC'
   ],
   [
     'CIN_low_vs_CIN_combined',
     'statistic'
   ],
   [
     'CIN_low_vs_CIN_combined',
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
    -rnk /home/chm2059/chm2059/0816_sam/results/RNAseq/DESeq2/${comparison}/${comparison}.${rnk}.rnk \
    -gmx /home/chm2059/chm2059/data/msigdb/${gmt} \
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
    set "${comparison}_${rnk}_${gmt}.xlsx" into get_output_out

  """
  python ${params.seqpy}/get_gsea_output.py \
    --indir ${fin} \
    --up_tab_name up_regulated \
    --down_tab_name down_regulated \
    --out ${comparison}_${rnk}_${gmt}
  """
}
