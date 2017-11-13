gzip -cd /Users/charlesmurphy/Desktop/Research/data/KEGG/170831/pathway.gz | sed -n "/PATHWAY_MAP hsa04010  MAPK signaling pathway/,/\/\/\//p" > MAPK.signlaing.pathway.txt &
gzip -cd /Users/charlesmurphy/Desktop/Research/data/KEGG/170831/pathway.gz | sed -n "/PATHWAY_MAP hsa04151  PI3K-Akt signaling pathway/,/\/\/\//p" > PI3K-Akt.signaling.pathway.txt &
