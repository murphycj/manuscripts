
java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/WT_vs_Rspondin_human/WT_vs_Rspondin_logPvalSignFC.rnk \
 -gmx /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/COAD/COAD_normal_vs_tumor/COAD_gene_signatures.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label WT_vs_Rspondin_TCGACOAD_logPvalSignFC \
 -include_only_symbols true \
 -make_sets true \
 -plot_top_x 40 \
 -rnd_seed timestamp \
 -set_max 5000 \
 -set_min 15 \
 -zip_report true \
 -out . \
 -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/WT_vs_Rspondin_human/WT_vs_Rspondin_statistic.rnk \
 -gmx /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/COAD/COAD_normal_vs_tumor/COAD_gene_signatures.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label WT_vs_Rspondin_TCGACOAD_statistic \
 -include_only_symbols true \
 -make_sets true \
 -plot_top_x 40 \
 -rnd_seed timestamp \
 -set_max 5000 \
 -set_min 15 \
 -zip_report true \
 -out . \
 -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/WT_vs_Rspondin_human/WT_vs_Rspondin_log2FC.rnk \
 -gmx /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/COAD/COAD_normal_vs_tumor/COAD_gene_signatures.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label WT_vs_Rspondin_TCGACOAD_log2FC \
 -include_only_symbols true \
 -make_sets true \
 -plot_top_x 40 \
 -rnd_seed timestamp \
 -set_max 5000 \
 -set_min 15 \
 -zip_report true \
 -out . \
 -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.1.jar -Xmx2048m xtools.gsea.GseaPreranked \
  -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/WT_vs_Rspondin_human/WT_vs_Rspondin_logPvalSignFC.rnk \
  -gmx ../../COAD.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label WT_vs_Rspondin_msigdbCOAD_logPvalSignFC \
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
 -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/WT_vs_Rspondin_human/WT_vs_Rspondin_statistic.rnk \
 -gmx ../../COAD.gmt \
 -collapse false \
 -mode Max_probe \
 -norm meandiv \
 -nperm 1000 \
 -scoring_scheme weighted \
 -rpt_label WT_vs_Rspondin_msigdbCOAD_statistic \
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
  -rnk /Users/charlesmurphy/Desktop/Research/dow_2016_5/results/RNAseq/DESeq2/WT_vs_Rspondin_human/WT_vs_Rspondin_log2FC.rnk \
  -gmx ../../COAD.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label WT_vs_Rspondin_msigdbCOAD_log2FC \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 40 \
  -rnd_seed timestamp \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out . \
  -gui false
