#Rscript APCMut_vs_APCWT.R
python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py --infile APCWT_vs_APCMUT_results.csv --reverseSign --ranking statistic --out APCWT_vs_APCMUT_results_statistic.rnk
python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py --infile APCWT_vs_APCMUT_results.csv --reverseSign --ranking logPvalSignFC --out APCWT_vs_APCMUT_results_logPvalSignFC.rnk
python ~/Desktop/Research/lib/seqpy/bin/deseq2rnk.py --infile APCWT_vs_APCMUT_results.csv --reverseSign --ranking log2FC --out APCWT_vs_APCMUT_results_log2FC.rnk
