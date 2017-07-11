
java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/APCWT_vs_APCMUT_human/APCWT_vs_APCMUT_results_logPvalSignFC.rnk \
 -gmx /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/COAD/COAD_normal_vs_tumor/COAD_gene_signatures.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label APCWT_vs_APCMUT_TCGACOAD_logPvalSignFC \
 -include_only_symbols true \
 -make_sets true \
 -plot_top_x 40 \
 -rnd_seed timestamp \
 -set_max 1500 \
 -set_min 15 \
 -zip_report true \
 -out . \
 -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/APCWT_vs_APCMUT_human/APCWT_vs_APCMUT_results_statistic.rnk \
 -gmx /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/COAD/COAD_normal_vs_tumor/COAD_gene_signatures.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label APCWT_vs_APCMUT_TCGACOAD_statistic \
 -include_only_symbols true \
 -make_sets true \
 -plot_top_x 40 \
 -rnd_seed timestamp \
 -set_max 1500 \
 -set_min 15 \
 -zip_report true \
 -out . \
 -gui false


java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/APCWT_vs_APCMUT_human/APCWT_vs_APCMUT_results_log2FC.rnk \
 -gmx /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/COAD/COAD_normal_vs_tumor/COAD_gene_signatures.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label APCWT_vs_APCMUT_TCGACOAD_log2FC \
 -include_only_symbols true \
 -make_sets true \
 -plot_top_x 40 \
 -rnd_seed timestamp \
 -set_max 1500 \
 -set_min 15 \
 -zip_report true \
 -out . \
 -gui false


java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
  -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/APCWT_vs_APCMUT_human/APCWT_vs_APCMUT_results_logPvalSignFC.rnk \
  -gmx ../../COAD.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label APCWT_vs_APCMUT_msigdbCOAD_logPvalSignFC \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 40 \
  -rnd_seed timestamp \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out . \
  -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/APCWT_vs_APCMUT_human/APCWT_vs_APCMUT_results_statistic.rnk \
 -gmx ../../COAD.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label APCWT_vs_APCMUT_msigdbCOAD_statistic \
 -include_only_symbols true \
 -make_sets true \
 -plot_top_x 40 \
 -rnd_seed timestamp \
 -set_max 500 \
 -set_min 15 \
 -zip_report true \
 -out . \
 -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
  -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/APCWT_vs_APCMUT_human/APCWT_vs_APCMUT_results_log2FC.rnk \
  -gmx ../../COAD.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label APCWT_vs_APCMUT_msigdbCOAD_log2FC \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 40 \
  -rnd_seed timestamp \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out . \
  -gui false
