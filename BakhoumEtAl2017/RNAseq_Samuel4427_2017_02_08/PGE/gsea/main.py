import pandas

genes = pandas.read_csv('../ensembl89.txt',sep='\t',header=None)
genes = genes[genes[1]=='1']
genes.index = genes[0]

deseq = pandas.read_csv('/Users/charlesmurphy/Desktop/Research/0816_sam/results/RNAseq_Samuel4427_2017_02_08/DESeq2_ensembl/primary_vs_metastasis/primary_vs_metastasis_results.csv',index_col=0)

common = list(set(deseq.index.tolist()).intersection(set(genes.index.tolist())))
deseq = deseq.loc[common,]

rnk = deseq['log2FoldChange']
rnk.index = genes.loc[common,4]
rnk.to_csv('primary_vs_metastasis.chr1.rnk',sep='\t')
deseq.index = genes.loc[deseq.index,4]

deseq.to_csv("de_genes.csv")
