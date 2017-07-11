#Rscript APCMut_vs_APCWT.R
python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py --infile COAD_normal_vs_tumor_results.csv --reverseSign --ranking statistic --out COAD_normal_vs_tumor_results_statistic.rnk
python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py --infile COAD_normal_vs_tumor_results.csv --reverseSign --ranking logPvalSignFC --out COAD_normal_vs_tumor_results_logPvalSignFC.rnk
python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py --infile COAD_normal_vs_tumor_results.csv --reverseSign --ranking log2FC --out COAD_normal_vs_tumor_results_log2FC.rnk
