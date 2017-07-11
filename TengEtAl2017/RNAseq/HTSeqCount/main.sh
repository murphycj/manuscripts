Rscript ~/Desktop/Research/lib/seqpy/bin/biomart_convert_gene_identifiers.R \
  -data HTSeq.gene.counts.csv \
  -d csv \
  -intype ensembl_gene_id \
  -outtype mgi_symbol \
  -dataset mmusculus_gene_ensembl \
  -out HTSeq.geneSymbol.counts.csv

python ~/Desktop/Research/lib/seqpy/bin/mouseSymbol2human.py --infile HTSeq.geneSymbol.vst.csv --out HTSeq.humanSymbols.vst.csv
