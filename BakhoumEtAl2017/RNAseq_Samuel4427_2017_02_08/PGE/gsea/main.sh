
java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.3.jar -Xmx2048m xtools.gsea.GseaPreranked \
  -rnk primary_vs_metastasis.chr1.rnk \
  -gmx /Users/charlesmurphy/Desktop/Research/data/msigdb/msigdb_v6.0_GMTs/h.all.v6.0.symbols.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label primary_vs_metastasis.h.all \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 80 \
  -rnd_seed timestamp \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out . \
  -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.3.jar -Xmx2048m xtools.gsea.GseaPreranked \
  -rnk primary_vs_metastasis.chr1.rnk \
  -gmx /Users/charlesmurphy/Desktop/Research/data/msigdb/msigdb_v6.0_GMTs/c2.cp.v6.0.symbols.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label primary_vs_metastasis.c2.cp \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 80 \
  -rnd_seed timestamp \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out . \
  -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.3.jar -Xmx2048m xtools.gsea.GseaPreranked \
  -rnk primary_vs_metastasis.chr1.rnk \
  -gmx /Users/charlesmurphy/Desktop/Research/data/msigdb/msigdb_v6.0_GMTs/c2.all.v6.0.symbols.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label primary_vs_metastasis.c2.all \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 80 \
  -rnd_seed timestamp \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out . \
  -gui false

java -cp /Users/charlesmurphy/Desktop/tools/gsea2-2.2.3.jar -Xmx2048m xtools.gsea.GseaPreranked \
  -rnk primary_vs_metastasis.chr1.rnk \
  -gmx /Users/charlesmurphy/Desktop/Research/data/msigdb/msigdb_v6.0_GMTs/c3.tft.v6.0.symbols.gmt \
  -collapse false \
  -mode Max_probe \
  -norm meandiv \
  -nperm 1000 \
  -scoring_scheme weighted \
  -rpt_label primary_vs_metastasis.c3.tft \
  -include_only_symbols true \
  -make_sets true \
  -plot_top_x 80 \
  -rnd_seed timestamp \
  -set_max 500 \
  -set_min 15 \
  -zip_report true \
  -out . \
  -gui false


python /Users/charlesmurphy/Desktop/Research/lib/seqpy/bin/get_gsea_output.py \
  --indir primary_vs_metastasis.c2.all.GseaPreranked.1502849305901/ \
  --up_tab_name up_regulated \
  --down_tab_name down_regulated \
  --leadingEdgeGenes \
  --allGenes \
  --leadingEdge \
  --out primary_vs_metastasis.c2.all
python /Users/charlesmurphy/Desktop/Research/lib/seqpy/bin/get_gsea_output.py \
  --indir primary_vs_metastasis.c2.cp.GseaPreranked.1502849272076/ \
  --up_tab_name up_regulated \
  --down_tab_name down_regulated \
  --leadingEdgeGenes \
  --allGenes \
  --leadingEdge \
  --out primary_vs_metastasis.c2.cp
python /Users/charlesmurphy/Desktop/Research/lib/seqpy/bin/get_gsea_output.py \
  --indir primary_vs_metastasis.h.all.GseaPreranked.1502849189599 \
  --up_tab_name up_regulated \
  --down_tab_name down_regulated \
  --leadingEdgeGenes \
  --allGenes \
  --leadingEdge \
  --out primary_vs_metastasis.h.all
python /Users/charlesmurphy/Desktop/Research/lib/seqpy/bin/get_gsea_output.py \
  --indir primary_vs_metastasis.c3.tft.GseaPreranked.1502849441582/ \
  --up_tab_name up_regulated \
  --down_tab_name down_regulated \
  --leadingEdgeGenes \
  --allGenes \
  --leadingEdge \
  --out primary_vs_metastasis.c3.tft
