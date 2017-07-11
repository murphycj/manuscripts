import pyensembl
import re
import os

data = pyensembl.EnsemblRelease(84,'mouse')

genes = []
genes += [i for i in data.gene_names() if re.search('^Igkv', i)]
genes += [i for i in data.gene_names() if re.search('^Igh', i)]
genes += [i for i in data.gene_names() if re.search('^Iglv', i)]
genes += [i for i in data.gene_names() if re.search('^H2-', i)]
genes += [i for i in data.gene_names() if re.search('^Klr', i)]

fout = open('excluded_genes.bed', 'w')
fout2 = open('excluded_genes.txt', 'w')
for i in genes:
    fout2.write(i + '\n')
    gg = data.genes_by_name(i)
    for gene in gg:
        fout.write(
            gene.contig + '\t' +
            str(gene.start) + '\t' +
            str(gene.end) + '\n'
        )
fout.close()
fout2.close()

os.system('~/Desktop/tools/bedtools2/bin/bedtools sort -i excluded_genes.bed > excluded_genes.sorted.bed ')
