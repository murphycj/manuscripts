import os
import pandas
import pyensembl

def parse_result_file(filename):
    mapping = []
    raw = []
    bed = []
    index=0


    data = open(filename,'r').read().split('\n')
    mapping.append(data[index].replace('<mapping>','').split(','))
    index+=1
    while data[index].find('</mapping><raw>')==-1:
        mapping.append(data[index].split(','))
        index+=1
    if data[index].find('</mapping><raw></raw>')==-1:
        raw.append(data[index].replace('</mapping><raw>','').split('\t'))
        index+=1
        while data[index].find('</raw><bed>')==-1:
            raw.append(data[index].split('\t'))
            index+=1
    else:
        index+=1

    while (index < len(data)) and (data[index].find('</bed>')==-1):
        bed.append(data[index].split('\t'))
        index+=1

    mapping = pandas.DataFrame(mapping,columns=['gene_id1','gene_id2','chrom','start','end','gene_id3'])
    raw = pandas.DataFrame(raw,columns=['chrom','start','end','pvalue','padj','num_differentially_expressed_genes','total_genes_in_region'])
    bed = pandas.DataFrame(bed,columns=['chrom','start','end','pval'])

    return mapping,raw,bed


def prep_data(comparisons,comparisons_params):
    for comparison in comparisons:
        outdir = comparison + '/pval' + str(comparisons_params[comparison]['pvalue']) + '_padj' + str(comparisons_params[comparison]['padj']) + '_log2FC' + str(comparisons_params[comparison]['log2'])
        if not os.path.exists(outdir):
            os.mkdir(outdir)
        deseq = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0816_sam/results/RNAseq_Samuel4427_2017_02_08/DESeq2_ensembl/" + comparison + "/" + comparison + "_results.csv")

        deseq2 = deseq[(deseq['padj']<=comparisons_params[comparison]['padj']) & (deseq['log2FoldChange'] >= comparisons_params[comparison]['log2']) & (deseq['pvalue']<=comparisons_params[comparison]['pvalue'])]
        deseq2['Unnamed: 0'].to_csv(outdir + '/' + comparison + '_up.txt',index=False,sep='\t')

        deseq2 = deseq[(deseq['padj']<=comparisons_params[comparison]['padj']) & (deseq['log2FoldChange']<= -comparisons_params[comparison]['log2']) & (deseq['pvalue']<=comparisons_params[comparison]['pvalue'])]
        deseq2['Unnamed: 0'].to_csv(outdir + '/' + comparison + '_down.txt',index=False,sep='\t')

def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def get_overlapping_genes(raw,mapping,genes):
    overlapping_genes = []
    mapping['chrom']=mapping['chrom'].str.upper()
    for i in raw.index:
        region_start = int(raw.loc[i,'start'])
        region_end = int(raw.loc[i,'end'])
        chrom = raw.loc[i,'chrom']
        mapping_sub = mapping[mapping['chrom']==chrom]

        subset = []
        for j in mapping_sub.index:
            j_start = int(mapping_sub.loc[j,'start'])
            j_end = int(mapping_sub.loc[j,'end'])
            if getOverlap([region_start,region_end],[j_start,j_end])>0:
                subset.append(True)
            else:
                subset.append(False)

        mapping_sub = mapping_sub[subset]['gene_id1'].tolist()
        mapping_sub = map(lambda x: x.upper(),mapping_sub)
        mapping_sub = ','.join(genes.loc[mapping_sub,'symbols'].tolist())
        overlapping_genes.append(mapping_sub)

    raw['differentially_expressed_genes'] = overlapping_genes
    raw['padj'] = raw['padj'].astype(float)
    raw = raw.sort_values(['padj'])
    return raw

def process_results(comparisons,comparisons_params):

    db = pyensembl.EnsemblRelease('75','human')

    for directory in comparisons:
        outdir = directory + '/pval' + str(comparisons_params[directory]['pvalue']) + '_padj' + str(comparisons_params[directory]['padj']) + '_log2FC' + str(comparisons_params[directory]['log2'])
        print outdir

        try:
            genes = pandas.read_table('./' + outdir + '/' + directory + '_up.txt',header=None)
        except:
            continue

        symbols = []
        for ii in genes[0].tolist():
            try:
                symbols.append(db.gene_by_id(ii).gene_name)
            except ValueError:
                symbols.append(ii)
        genes['symbols'] = symbols
        genes.index = genes[0]

        mapping,raw,bed = parse_result_file(
            './' + outdir + '/' + directory + '_results_up.txt'
        )
        raw = get_overlapping_genes(raw,mapping,genes)
        raw_up = raw[raw['total_genes_in_region'].astype(int)>3]
        raw_up['cna'] = ['amplification']*raw_up.shape[0]

        try:
            genes = pandas.read_table('./' + outdir + '/' + directory + '_down.txt',header=None)
        except:
            continue
        symbols = []
        for ii in genes[0].tolist():
            try:
                symbols.append(db.gene_by_id(ii).gene_name)
            except ValueError:
                symbols.append(ii)
        genes['symbols'] = symbols
        genes.index = genes[0]

        mapping,raw,bed = parse_result_file(
            './' + outdir + '/' + directory + '_results_down.txt'
        )
        raw = get_overlapping_genes(raw,mapping,genes)
        raw_down = raw[raw['total_genes_in_region'].astype(int)>3]
        raw_down['cna'] = ['deletion']*raw_down.shape[0]

        raw = raw_down.append(raw_up)
        raw = raw[['chrom','start','end','pvalue','padj','cna','num_differentially_expressed_genes','total_genes_in_region','differentially_expressed_genes']]
        raw = raw.sort_values('padj')

        raw.to_excel('./' + outdir + '/' + directory + '_results.xlsx',index=False)

def run_data(fout,comparisons,comparisons_params):
    for comparison in comparisons:
        outdir = comparison + '/pval' + str(comparisons_params[comparison]['pvalue']) + '_padj' + str(comparisons_params[comparison]['padj']) + '_log2FC' + str(comparisons_params[comparison]['log2'])
        fout.write(
            'perl pge.pl -q ./%s/%s_up.txt -r ensembl89 > ./%s/%s_results_up.txt\n' %
            (
                outdir,
                comparison,
                outdir,
                comparison
            )
        )
        fout.write(
            'perl pge.pl -q ./%s/%s_down.txt -r ensembl89 > ./%s/%s_results_down.txt\n' %
            (
                outdir,
                comparison,
                outdir,
                comparison
            )
        )

comparisons = [
    'primary_vs_metastasis','BoneMet_vs_BrainMet',
    'CINhighPrimary_vs_CINhighMetastasis','CINLowMets_vs_CINHighMets',
    'CINLowPrimary_vs_CINMedHighPrimary','Control_vs_STINGshRNA',
    'CINlow_vs_CINmedhigh']
comparisons=['1','2','3','23','24','31','32','45']


fout = open('pge_scripts.sh','w')
for pv in [0.05]:
    #for log2 in [0,0.25,0.5]:
    for log2 in [0]:
        #for padj in [0.05,0.1,0.25,1]:
        for padj in [0.05]:
            comparisons_params = {}
            comparisons_params['1'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['2'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['3'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['23'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['24'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['31'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['32'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['45'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['primary_vs_metastasis'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['BoneMet_vs_BrainMet'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINhighPrimary_vs_CINhighMetastasis'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINLowMets_vs_CINHighMets'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINLowPrimary_vs_CINMedHighPrimary'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['Control_vs_STINGshRNA'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINlow_vs_CINmedhigh'] = {'padj':padj,'pvalue':pv,'log2':log2}

            #prep_data(comparisons,comparisons_params)
            #run_data(fout,comparisons,comparisons_params)
fout.close()
#exit()
for pv in [0.05]:
    #for log2 in [0,0.25,0.5]:
    for log2 in [0]:
        #for padj in [0.05,0.1,0.25,1]:
        for padj in [0.05]:
            comparisons_params = {}
            comparisons_params['1'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['2'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['3'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['23'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['24'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['31'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['32'] = {'padj':padj,'pvalue':pv,'log2':log2}
            comparisons_params['45'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['primary_vs_metastasis'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['BoneMet_vs_BrainMet'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINhighPrimary_vs_CINhighMetastasis'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINLowMets_vs_CINHighMets'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINLowPrimary_vs_CINMedHighPrimary'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['Control_vs_STINGshRNA'] = {'padj':padj,'pvalue':pv,'log2':log2}
            #comparisons_params['CINlow_vs_CINmedhigh'] = {'padj':padj,'pvalue':pv,'log2':log2}

            process_results(comparisons,comparisons_params)
