Rscript ~/Desktop/Research/lib/seqpy/bin/biomart_convert_gene_identifiers.R \
  -data genes-fpkm.csv \
  -d csv \
  -intype ensembl_gene_id \
  -outtype mgi_symbol \
  -dataset mmusculus_gene_ensembl \
  -out gene-symbols-fpkm.csv
python ~/Desktop/Research/lib/seqpy/bin/mouseSymbol2human.py \
  --infile gene-symbols-fpkm.csv \
  --out gene-humanSymbols-fpkm.csv
Rscript ~/Desktop/Research/lib/seqpy/bin/biomart_convert_gene_identifiers.R \
  -data genes-humanSymbols-fpkm.csv \
  -d csv \
  -intype hgnc_symbol \
  -outtype entrezgene \
  -dataset hsapiens_gene_ensembl \
  -out genes-humanEntrez-fpkm.csv
