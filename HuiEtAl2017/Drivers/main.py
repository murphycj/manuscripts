import re
import pandas
import pyensembl
from Bio import AlignIO
from Bio.SeqUtils import seq1
import os
import pickle

def align(aname,aseq,bname,bseq):
    afile = './tmp/' + aname + '.fa'
    bfile = './tmp/' + bname + '.fa'
    outfile = './tmp/' + aname + '-' + bname + '.needle'

    fout = open(afile,'w')
    fout.write('>' + aname + '\n')
    fout.write(aseq)
    fout.close()

    fout = open(bfile,'w')
    fout.write('>' + bname + '\n')
    fout.write(bseq)
    fout.close()

    os.system(
        'needle -asequence %s -bsequence %s -datafile EBLOSUM62 -gapopen 10 -gapextend 0.5 -outfile %s' %
        (
            afile,
            bfile,
            outfile
        )
    )

    alignments = list(AlignIO.parse(outfile, "emboss"))

    os.system('rm ' + afile)
    os.system('rm ' + bfile)
    os.system('rm ' + outfile)

    if len(alignments) != 1:
        return None
    else:
        return alignments[0]

def create_align_db():

    alignments = dict()
    human_db = pyensembl.EnsemblRelease(86,'human')
    mouse_db = pyensembl.EnsemblRelease(86,'mouse')

    homology = pandas.read_csv("mart_export.txt")
    homology = homology[~pandas.isnull(homology['Mouse gene stable ID'])]
    homologs = dict(zip(homology['Mouse gene stable ID'],homology['Gene stable ID']))

    mutations = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/Mutations/Varscan-paired/somatic.nodbsnp.varscan.filtered.primaryControl.csv")
    for i in mutations.index:

        #get the canonical mouse transcript
        mouse_transcript = mutations.loc[i,'TRANSCRIPT'].split('.')[0]
        mouse_transcript = mouse_db.transcript_by_id(mouse_transcript)

        # determine if there is a human homolog

        if mouse_transcript.gene_id in homologs:
            human_gene = human_db.gene_by_id(homologs[mouse_transcript.gene_id])

            #get longest CDS transcript

            max_length = 0
            human_canonical_transcript = None

            for human_transcript in human_gene.transcripts:
                if human_transcript.protein_sequence is not None and len(human_transcript.protein_sequence)>max_length:
                    human_canonical_transcript=human_transcript
                    max_length = len(human_transcript.protein_sequence)

            # just skip ENST00000589042 because it is really long, and I know
            # none of mouse tumors have mutations in it

            if human_canonical_transcript is None or human_canonical_transcript.id == 'ENST00000589042':
                continue

            # align the protein sequences

            if human_canonical_transcript.id + '-' + mouse_transcript.id not in alignments:
                print human_canonical_transcript.id + '-' + mouse_transcript.id, len(human_canonical_transcript.protein_sequence), len(mouse_transcript.protein_sequence)
                alignments[human_canonical_transcript.id + '-' + mouse_transcript.id] = align(
                    human_canonical_transcript.id,
                    human_canonical_transcript.protein_sequence,
                    mouse_transcript.id,
                    mouse_transcript.protein_sequence
                )

    pickle.dump(alignments,open('alignments.p','wb'))

def parse_pathway(filename):
    genes = []
    parsing_genes = False
    for line in open(filename,'r').readlines():
        if line.find('GENE')!=-1:
            parsing_genes=True
        elif line.find('COMPOUND')!=-1:
            parsing_genes=False

        if parsing_genes:
            gene = re.findall('(?<=[0-9]  )(.*)(?=;)',line)
            if len(gene)!=1:
                print 'error %s' % gene
            genes.append(gene[0])

    return genes


def fetch_kegg_pathway_genes():

    genes = {'YAP1':['PI3K']}

    for gg in parse_pathway('/Users/charlesmurphy/Desktop/Research/0914_hui/results/Drivers/MAPK.signlaing.pathway.txt'):
        if gg not in genes:
            genes[gg] = ['MAPK']
        else:
            print '%s is already there' % gg
    for gg in parse_pathway('/Users/charlesmurphy/Desktop/Research/0914_hui/results/Drivers/PI3K-Akt.signaling.pathway.txt'):
        if gg not in genes:
            genes[gg] = ['PI3K']
        else:
            genes[gg].append('PI3K')


    genes2 = genes.keys()
    genes2.sort()

    fout = open('pathway_genes.csv','w')
    fout.write('gene,pathway\n')
    for gg in genes2:
        fout.write('%s,%s\n' % (gg,';'.join(genes[gg])))
    fout.close()


def convert_mutations():
    chroms = [str(i) for i in range(1,20)] + ['X','Y']

    # human and mouse homologs

    homology = pandas.read_csv("mart_export.txt")
    homology = homology[~pandas.isnull(homology['Mouse gene stable ID'])]
    homologs = dict(zip(homology['Mouse gene stable ID'],homology['Gene stable ID']))

    human_db = pyensembl.EnsemblRelease(86,'human')
    mouse_db = pyensembl.EnsemblRelease(86,'mouse')

    # alignments

    alignments = pickle.load(open('alignments.p','rb'))
    alignment_maps = {}
    for i in alignments.keys():
        alignment_maps[i.split('-')[1]] = i

    mutations = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/Mutations/Varscan-paired/somatic.nodbsnp.varscan.filtered.primaryControl.csv")
    hotspots = pandas.read_table('/Users/charlesmurphy/Desktop/Research/data/papers/hotspotMutations/060717.hotspots.txt',sep='\t')
    ishotspot = []
    subsitution = []
    human_genes = []

    for i in mutations.index:

        #get the mouse transcript and determine if it had an alignment

        transcript = mutations.loc[i,'TRANSCRIPT'].split('.')[0]
        if transcript not in alignment_maps:
            ishotspot.append('')
            subsitution.append('')
            human_genes.append('')
            continue

        human_gene = human_db.transcript_by_id(alignment_maps[transcript].split('-')[0]).gene.name
        human_genes.append(human_gene)

        AA = mutations.loc[i,'AMINO ACID CHANGE']
        bp_change = mutations.loc[i,'BASE PAIR CHANGE']
        effect = mutations.loc[i,'EFFECT']

        if pandas.isnull(AA) and effect.find('splice')!=-1:

            #only continute if it is an actualy amino acid change
            #splice site variants don't have any entry for AA, but I manually
            #checked the list in OncoKB and saw no splice site variants

            ishotspot.append('')
            subsitution.append(bp_change + ' splice variant')
            continue
        elif pandas.isnull(AA) and effect.find('5_prime_UTR_premature_start_codon_gain_variant')!=-1:
            ishotspot.append('')
            subsitution.append(bp_change + ' 5_prime_UTR_premature_start_codon_gain_variant')
            continue
        elif AA.find('*')!=-1 or AA.find('fs')!=-1 or AA=='p.Met1?':
            ishotspot.append('')
            subsitution.append(AA + ' truncating')
            continue
        elif pandas.isnull(AA):
            ishotspot.append('')
            subsitution.append('')
            print 'not sure what variant'
            import pdb; pdb.set_trace()
            continue

        # get the amino acid position and letter change

        AA_pos = int(re.findall('[0-9]+',AA)[0])
        AA = AA.replace('p.','')
        AA = AA.replace(str(AA_pos),'')
        AA = seq1(AA)

        # determine if there was an alignment

        alignment = alignments[alignment_maps[transcript]]

        if alignment is None:
            ishotspot.append('')
            subsitution.append(AA[0] + str(human_pos) + AA[1])
            continue

        # if there was an alignment, get the amino acid position in the human gene

        mouse_pos = 0
        mouse_index = 0
        mouse_aa = ''
        for aa in alignment[1]:
            mouse_index += 1
            if aa=='X':
                continue
            if aa != '-':
                mouse_pos += 1
            if mouse_pos==AA_pos:
                mouse_aa = alignment[1][mouse_index-1]
                break

        if mouse_aa != AA[0]:
            print 'boo'
            import pdb; pdb.set_trace()

        human_pos = 0
        human_index = 0
        human_aa = ''
        for aa in alignment[0]:
            human_index += 1
            if aa=='X':
                continue
            if aa != '-':
                human_pos += 1
            if human_index==mouse_index:
                human_aa = alignment[0][human_index-1]
                break

        subsitution.append(AA[0] + str(human_pos) + AA[1])

        # determine if it is in a hotspot

        hotspots_tmp = hotspots[hotspots['Gene']==human_gene]
        if hotspots_tmp.shape[0]>0:
            positions = map(lambda x: int(re.findall('[0-9]+',x)[0]), hotspots_tmp['Residue'].tolist())
            positions = [human_pos == i for i in positions]
            if sum(positions) > 0:
                hotspots_tmp = hotspots_tmp[positions]['Residue'].tolist()[0]
                ishotspot.append(human_gene + ' (' + hotspots_tmp + ')')
            else:
                ishotspot.append('')
        else:
            ishotspot.append('')

    mutations_per_sample = {}

    samples = mutations.columns.tolist()[10:]
    for sample in samples:
        mutations_per_sample[sample] = {}
        #mutations_per_sample[sample]['genes'] = []
        mutations_per_sample[sample]['truncating'] = []
        mutations_per_sample[sample]['premature'] = []
        mutations_per_sample[sample]['splice'] = []
        mutations_per_sample[sample]['missense'] = []
        mutations_per_sample[sample]['hotspot'] = []
    for i in range(0, mutations.shape[0]):
        for ss in samples:
            if mutations.loc[i,ss]!='-':
                #mutations_per_sample[ss]['genes'].append(human_genes[i])
                if ishotspot[i]!='':
                    mutations_per_sample[ss]['hotspot'].append(ishotspot[i])
                if subsitution[i].find('truncating')!=-1:
                    mutations_per_sample[ss]['truncating'].append(human_genes[i] + ' (' + subsitution[i] + ')')
                elif subsitution[i].find('splice')!=-1:
                    mutations_per_sample[ss]['splice'].append(human_genes[i] + ' (' + subsitution[i] + ')')
                elif subsitution[i].find('5_prime_UTR_premature_start_codon_gain_variant')!=-1:
                    mutations_per_sample[ss]['premature'].append(human_genes[i] + ' (' + subsitution[i] + ')')
                elif subsitution[i]!='':
                    mutations_per_sample[ss]['missense'].append(human_genes[i] + '-' + subsitution[i])

    return mutations_per_sample

def annotate_cnas():

    cnas_per_sample = {}

    cnas = pandas.read_csv('/Users/charlesmurphy/Desktop/Research/0914_hui/results/WEX/CNV/CNVkit/gene-sample-cn-humanGeneSymbols.csv',index_col=0)

    for sample in cnas.columns:
        cnas_per_sample[sample] = {}
        cnas_per_sample[sample]['AMP'] = cnas[cnas[sample]>3].index.tolist()
        cnas_per_sample[sample]['DEL'] = cnas[cnas[sample]==0].index.tolist()

    return cnas_per_sample

def annotate_expression(samples,threshold):

    expression_per_sample = {}
    fpkm = pandas.read_csv("/Users/charlesmurphy/Desktop/Research/0914_hui/results/RNAseq/Cufflinks/genes-humanSymbols-fpkm.csv",index_col=0)
    fpkm[fpkm<=0.1]=0
    fpkm=fpkm.fillna(0)

    common = set(samples['sample_ID'].tolist()).intersection(set(fpkm.columns.tolist()))
    fpkm = fpkm[list(common)]

    for i in fpkm.index:
        fpkm.loc[i,:] = (fpkm.loc[i,:]-fpkm.loc[i,:].mean())/fpkm.loc[i,:].std()

    for sample in common:
        expression_per_sample[sample] = {}
        expression_per_sample[sample]['UP'] = fpkm[fpkm[sample]>=threshold].index.tolist()
        expression_per_sample[sample]['DOWN'] = fpkm[fpkm[sample]<=-threshold].index.tolist()

    return expression_per_sample

def annotate_fusions(samples):

    fusions_per_sample = {}

    homology = pandas.read_csv("mart_export.txt")
    homology = homology[~pandas.isnull(homology['Mouse gene stable ID'])]
    homologs = dict(zip(homology['Mouse gene stable ID'],homology['Gene stable ID']))

    fusions = pandas.read_csv('/Users/charlesmurphy/Desktop/Research/0914_hui/results/Fusions/primaryControl/fusions-noNormal.noReadthrough.n2.csv')

    for sample in samples['sample_ID'].tolist():
        tmp = fusions[fusions['sample']==sample]
        fusions_per_sample[sample] = {}
        fusions_per_sample[sample]['genes'] = []
        fusions_per_sample[sample]['fusions'] = []
        if tmp.shape[0]>0:
            for i in tmp.index:
                gene1 = tmp.loc[i,'Gene_1_id.5end_fusion_partner.']
                gene2 = tmp.loc[i,'Gene_2_id.3end_fusion_partner.']
                gene1_human = gene2_human = ''

                if gene1 in homologs:
                    tmp2 = homology[homology['Gene stable ID']==homologs[gene1]]
                    gene1_human = tmp2['Gene name'].tolist()[0]
                    if gene1_human=='':
                        import pdb; pdb.set_trace()
                    fusions_per_sample[sample]['genes'].append(gene1_human)

                if gene2 in homologs:
                    tmp2 = homology[homology['Gene stable ID']==homologs[gene2]]
                    gene2_human = tmp2['Gene name'].tolist()[0]
                    if gene2_human=='':
                        import pdb; pdb.set_trace()
                    fusions_per_sample[sample]['genes'].append(gene2_human)

                if gene1_human != '' and gene2_human != '':
                    fusions_per_sample[sample]['fusions'].append(gene1_human + ';' + gene2_human)

    for sample in fusions_per_sample.keys():
        fusions_per_sample[sample]['fusions'] = list(set(fusions_per_sample[sample]['fusions']))
        fusions_per_sample[sample]['genes'] = list(set(fusions_per_sample[sample]['genes']))

    return fusions_per_sample

def prep_oncokb_data(filename):

    results = {}

    oncokb = pandas.read_table(filename,sep='\t')
    oncokb = oncokb[(oncokb['Oncogenicity']!='Inconclusive') & (oncokb['Oncogenicity']!='Likely Neutral') & (oncokb['Oncogenicity']!='null')]
    results['overexpressed'] = set(oncokb[oncokb['Alteration']=='Overexpression']['Gene'].tolist())
    results['fusion'] = set(oncokb[oncokb['Alteration']=='Fusions']['Gene'].tolist())
    oncokb = oncokb[(oncokb['Alteration']!='Overexpression') & (oncokb['Alteration']!='Fusions')]
    results['fusion-pair'] = set(oncokb[oncokb['Alteration'].str.contains('Fusion')]['Alteration'].str.replace(' Fusion','').tolist())
    oncokb = oncokb[~oncokb['Alteration'].str.contains('Fusion')]
    results['AMP'] = set(oncokb[oncokb['Alteration']=='Amplification']['Gene'].tolist())
    results['DEL'] = set(oncokb[oncokb['Alteration']=='Deletion']['Gene'].tolist())
    results['truncating'] = set(oncokb[oncokb['Alteration']=='Truncating Mutations']['Gene'].tolist())

    # remove alteration types not measured

    remove = [
        'Promoter Mutations','Promoter Mutations',
        'Epigenetic Silencing','3\' Deletion','Hypermethylation',
        'Amplification','Truncating Mutations','Deletion',
        'Promoter Hypermethylation'
    ]
    oncokb = oncokb[~oncokb['Alteration'].isin(remove)]

    results['mutations'] = {}
    for g in oncokb.groupby(['Gene']):
        results['mutations'][g[0]] = []
        for i in g[1]['Alteration']:
            results['mutations'][g[0]].append(i)

    return results

def annotate_with_oncoKB():

    pi3k_mapk_genes = pandas.read_csv('pathway_genes.csv')
    pi3k_mapk_genes = dict(zip(pi3k_mapk_genes['gene'].tolist(),pi3k_mapk_genes['pathway'].tolist()))

    samples = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
    samples = samples[(samples['byl']==0) & (samples['bkm']==0) & (samples['bmn']==0)]
    samples = samples[(samples['tissue']!='liver') & (samples['tissue']!='normal breast') & (samples['tissue']!='tail')]
    samples['sample_name'] = samples['sample_name'].astype(str)
    samples.index = samples['sample_ID']
    samples_rnaseq = samples[samples['sequencing_type']=='RNAseq']
    samples_wex = samples[samples['sequencing_type']=='WEX']

    mutations = convert_mutations()
    fusions = annotate_fusions(samples_rnaseq)
    cnas = annotate_cnas()
    expression = annotate_expression(samples_rnaseq, 1.5)

    data = {
        i:{
            'Hotspot(s)':[],
            'Oncogenic':[],
            'Drug(s)':[],
            'Actionable':[],
            'PI3K':[],
            'MAPK':[]
        } for i in samples['sample_ID'].tolist()
    }

    oncokb = prep_oncokb_data("/Users/charlesmurphy/Desktop/Research/data/OncoKB/170831-allAnnotatedVariants.txt")
    oncokb_actionable_table = pandas.read_table("/Users/charlesmurphy/Desktop/Research/data/OncoKB/170831-allActionableVariants.txt",sep='\t')

    # manually check for mutations

    data['HL107']['Oncogenic'] = ['JAK1-L910P']
    data['HL107']['PI3K'] = ['JAK1 (L910P)']
    data['HL153']['Oncogenic'] = ['JAK1-L910P']
    data['HL153']['PI3K'] = ['JAK1 (L910P)']
    data['HL160']['Oncogenic'] = ['JAK1-F958C']
    data['HL160']['PI3K'] = ['JAK1 (F958C)']
    data['HL23']['Oncogenic'] = ['JAK1-F958C']
    data['HL23']['PI3K'] = ['JAK1 (F958C)']
    data['HL-WES-66']['Oncogenic'] = ['HRAS-Q61K']
    data['HL-WES-66']['PI3K'] = ['HRAS (Q61K)']
    data['HL-WES-66']['MAPK'] = ['HRAS (Q61K)']
    data['HL162']['Oncogenic'] = ['HRAS-Q61K']
    data['HL162']['PI3K'] = ['HRAS (Q61K)']
    data['HL162']['MAPK'] = ['HRAS (Q61K)']
    data['HL166']['Oncogenic'] = ['KRAS-Q61H']
    data['HL166']['PI3K'] = ['KRAS (Q61H)']
    data['HL166']['MAPK'] = ['KRAS (Q61H)']
    data['HL166']['Actionable'] = ['KRAS-Q61H']
    data['HL166']['Drug(s)'] = [
        'Docetaxel + Trametinib','Palbociclib + Trametinib','Erlotinib + Binimetinib',
        'Selumetinib','Abemaciclib + Trametinib','Binimetinib','Ribociclib + Trametinib',
        'Trametinib','LY3214996','KO-947','GDC-0994','Palbociclib + Trametinib',
        'Binimetinib + Ribociclib','Cobimetinib + Atezolizumab'
    ]

    additional_genes = ['FGFR2','MET','YAP1']

    for sample in samples['sample_ID'].tolist():

        if sample in expression:

            for gene in expression[sample]['UP']:
                if gene in oncokb['overexpressed'] or gene in additional_genes:
                    data[sample]['Oncogenic'].append(gene + ' (overexpressed)')

                if gene in additional_genes:
                    for i in pi3k_mapk_genes[gene].split(';'):
                        data[sample][i].append(gene + ' (overexpressed)')

        # mutations

        if sample in mutations:

            for gene in mutations[sample]['truncating']:
                if gene.split(' ')[0] in oncokb['truncating']:

                    tmp = oncokb_actionable_table[(oncokb_actionable_table['Gene']==gene.split(' ')[0]) & (oncokb_actionable_table['Alteration']=='Truncating Mutations')]
                    if tmp.shape[0]>0:
                        data[sample]['Actionable'].append(gene)
                        data[sample]['Drug(s)'] += list(set((', '.join(tmp['Drugs(s)'].tolist())).split(', ')))


                    data[sample]['Oncogenic'].append(gene)

                    tmp = oncokb_actionable_table[(oncokb_actionable_table['Gene']==gene.split(' ')[0]) & (oncokb_actionable_table['Alteration']=='Oncogenic Mutations')]
                    if tmp.shape[0]>0:
                        data[sample]['Actionable'].append(gene)
                        data[sample]['Drug(s)'] += list(set((', '.join(tmp['Drugs(s)'].tolist())).split(', ')))

                if gene.split(' ')[0] in pi3k_mapk_genes:
                    for i in pi3k_mapk_genes[gene.split(' ')[0]].split(';'):
                        data[sample][i].append(gene)

            for ii in ['splice','missense','premature']:
                for gene in mutations[sample][ii]:
                    if gene.split(' ')[0] in pi3k_mapk_genes:
                        for i in pi3k_mapk_genes[gene.split(' ')[0]].split(';'):
                            data[sample][i].append(gene)

        # CNAs

        if sample in cnas:
            for gene in cnas[sample]['AMP']:
                if gene in oncokb['AMP'] or gene in additional_genes:
                    data[sample]['Oncogenic'].append(gene + ' (AMP)')

                    tmp = oncokb_actionable_table[(oncokb_actionable_table['Gene']==gene) & (oncokb_actionable_table['Alteration']=='Amplification')]
                    if tmp.shape[0]>0:
                        data[sample]['Actionable'].append(gene + ' (AMP)')
                        data[sample]['Drug(s)'] += list(set((', '.join(tmp['Drugs(s)'].tolist())).split(', ')))

                    if gene in additional_genes:
                        data[sample]['Actionable'].append(gene + ' (AMP)')

                if gene in pi3k_mapk_genes:
                    for i in pi3k_mapk_genes[gene].split(';'):
                        data[sample][i].append(gene + ' (AMP)')

            for gene in cnas[sample]['DEL']:
                if gene in oncokb['DEL'] or gene in additional_genes:
                    data[sample]['Oncogenic'].append(gene + ' (DEL)')

                    tmp = oncokb_actionable_table[(oncokb_actionable_table['Gene']==gene) & (oncokb_actionable_table['Alteration']=='Deletion')]
                    if tmp.shape[0]>0:
                        data[sample]['Actionable'].append(gene + ' (DEL)')
                        data[sample]['Drug(s)'] += list(set((', '.join(tmp['Drugs(s)'].tolist())).split(', ')))

                    if gene in additional_genes:
                        data[sample]['Actionable'].append(gene + ' (DEL)')

                if gene in pi3k_mapk_genes:
                    for i in pi3k_mapk_genes[gene].split(';'):
                        data[sample][i].append(gene + ' (DEL)')

        # Fusions

        if sample in fusions:
            for gene in fusions[sample]['genes']:
                if gene in oncokb['fusion']:
                    data[sample]['Oncogenic'].append(gene + ' (fusion)')

                    tmp = oncokb_actionable_table[(oncokb_actionable_table['Gene']==gene) & (oncokb_actionable_table['Alteration']=='Fusions')]
                    if tmp.shape[0]>0:
                        data[sample]['Actionable'].append(gene + ' (fusion)')
                        data[sample]['Drug(s)'] += list(set((', '.join(tmp['Drugs(s)'].tolist())).split(', ')))

                if gene in pi3k_mapk_genes:
                    for i in pi3k_mapk_genes[gene].split(';'):
                        data[sample][i].append(gene + ' (fusion)')

            for gene in fusions[sample]['fusions']:
                if gene in oncokb['fusion-pair']:
                    data[sample]['Oncogenic'].append(gene + ' (fusion)')

                    tmp = oncokb_actionable_table[oncokb_actionable_table['Alteration']==gene + ' Fusion']
                    if tmp.shape[0]>0:
                        data[sample]['Actionable'].append(gene + ' (fusion)')
                        data[sample]['Drug(s)'] += list(set((', '.join(tmp['Drugs(s)'].tolist())).split(', ')))

                if gene in pi3k_mapk_genes:
                    for i in pi3k_mapk_genes[gene].split(';'):
                        data[sample][i].append(gene + ' (fusion)')

    fout = open('tumor-drivers.tsv','w')
    fout.write('Primary_tumor_number\tSamples\tBRCA1\tPI3K_alterations\tMAPK_alterations\tHotspots\tOncogenic_alterations\tActionable_alterations\tDrug(s)\n')
    for group in samples.groupby(['mouse','primary_tumor_index']):
        if str(group[0][0]) in ['517','1783','1867','2029']:
            primary_tumor_number = str(group[0][0]) + '-' + str(group[0][1])
        else:
            primary_tumor_number = str(group[0][0])
        tmp_samples = ', '.join((group[1]['sample_ID'] + ' (' + group[1]['sequencing_type'] + ')').tolist())

        PI3K_alterations = []
        MAPK_alterations = []
        Hotspots = []
        Oncogenic_alterations = []
        Actionable_alterations = []
        Drugs = []
        for ss in group[1]['sample_ID']:
            if ss in data:
                PI3K_alterations += data[ss]['PI3K']
                MAPK_alterations += data[ss]['MAPK']
                Oncogenic_alterations += data[ss]['Oncogenic']
                Actionable_alterations += data[ss]['Actionable']
                Drugs += data[ss]['Drug(s)']
            if ss in mutations:
                Hotspots += mutations[ss]['hotspot']


        fout.write(
            '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %
            (
                primary_tumor_number,
                tmp_samples,
                group[1]['brca1'].tolist()[0],
                ', '.join(set(PI3K_alterations)),
                ', '.join(set(MAPK_alterations)),
                ', '.join(set(Hotspots)),
                ', '.join(set(Oncogenic_alterations)),
                ', '.join(set(Actionable_alterations)),
                ', '.join(set(Drugs))
            )
        )
    fout.close()

#create_align_db()
#fetch_kegg_pathway_genes()
annotate_with_oncoKB()
