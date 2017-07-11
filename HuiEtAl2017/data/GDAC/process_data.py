import pandas
import mouse2human

data = pandas.read_table(
    "../BRCA.rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.data.txt",
    sep='\t'
)


genes = data['Hybridization REF'].tolist()[1:]
new_genes = []
indices_to_keep = []
for i in data.index:
    temp = data.ix[i]['Hybridization REF'].split('|')[0]
    if temp!='?' and temp not in new_genes:
        new_genes.append(temp)
        indices_to_keep.append(i)

data = data.ix[indices_to_keep]
data.index = new_genes

#print counts

counts_columns = data.ix['gene']=='raw_counts'
counts = data[data.columns[counts_columns]]

samples = counts.columns.tolist()
simple_sample_names = []
for i in samples:
    simple_sample_names.append('-'.join(i.split('-')[0:4]))
counts.columns = simple_sample_names
counts = counts.ix[counts.index.tolist()[1:]]
counts.to_csv("BRCA1.counts.csv")

#print RPKM

RPKM_columns = data.ix['gene']=='RPKM'
RPKM = data[data.columns[RPKM_columns]]

samples = RPKM.columns.tolist()
simple_sample_names = []
for i in samples:
    simple_sample_names.append('-'.join(i.split('-')[0:4]))
RPKM.columns = simple_sample_names

RPKM = RPKM.ix[RPKM.index.tolist()[1:]]
RPKM.to_csv("BRCA.rpkm.csv")
