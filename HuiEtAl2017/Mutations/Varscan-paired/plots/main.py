from collections import Counter
import vcf
import pandas
from seqpy.parsers import *

def count_by_gene():
    samples = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
    sample_map = dict(zip(samples['sample_ID'].tolist(),samples['mouse'].tolist()))

    primaries = list(set(samples['mouse'].tolist()))

    vcf_in = vcf.Reader(open('../somatic.nodbsnp.varscan.filtered.vcf','r'))

    genes={}
    temp = [
        'total_mutations','n_synonymous','n_missense','n_spice','n_frameshift',
        'n_truncating','n_inframe_indel','total_sample','total_primary'
    ]

    for variant in vcf_in:

        variant_info = SnpEffInfo(variant.INFO)
        called_samples = []
        for s in variant.samples:
            if s.called:
                called_samples.append(s.sample)

        n = len(called_samples)

        primary_count = {i:False for i in primaries}
        for s in called_samples:
            primary_count[sample_map[s]]=True

        n_primary = sum(map(lambda x: 1 if x else 0, primary_count.values()))

        if variant_info.has_ann():
            gene_effects = variant_info.get_effects_by_gene()

            for gene, effects in gene_effects.items():
                if gene not in genes:
                    genes[gene] = {'total':n,'total_primary':n_primary}
                else:
                    genes[gene]['total'] += n
                    genes[gene]['total_primary'] += n_primary

                for e in effects:
                    if e not in genes[gene]:
                        genes[gene][e] = n_primary
                    else:
                        genes[gene][e] += n_primary

    effects = EFFECTS.keys()
    effects.remove('sequence_feature + exon_loss_variant')
    effects.remove('3_prime_UTR_variant')
    effects.remove('miRNA')
    effects.remove('exon_variant')
    effects.remove('5_prime_UTR_variant')
    effects.remove('conserved_intergenic_variant')
    effects.remove('non_coding_exon_variant')
    effects.remove('intragenic_variant')
    effects.remove('downstream_gene_variant')
    effects.remove('transcript_variant')
    effects.remove('chromosome')
    effects.remove('coding_sequence_variant')
    effects.remove('gene_variant')
    effects.remove('upstream_gene_variant')
    effects.remove('rare_amino_acid_variant')
    effects.remove('intergenic_region')
    effects.remove('conserved_intron_variant')
    effects.remove('regulatory_region_variant')
    effects.remove('intron_variant')
    effects.remove('exon_loss_variant')
    effects.remove('3_prime_UTR_truncation + exon_loss')



    fout = open('gene.count.csv','w')
    fout.write('gene,total_mutations,total_primary,' + ','.join(effects) + '\n')
    for gene,stats in genes.items():
        fout.write(gene + ',' + str(stats['total']) + ',' + str(stats['total_primary']))

        for i in effects:
            if i in stats:
                fout.write(',' + str(stats[i]))
            else:
                fout.write(',0')
        fout.write('\n')
    fout.close()

def compute_shared_within_primary():
    samples = pandas.read_excel('/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx','samples')
    sample_map = dict(zip(samples['sample_ID'].tolist(),samples['mouse'].tolist()))

    ps = [1720,1513,1259,1397,1512,1367,1616,1460,1221,1536,1461,1413,1795,1660,1661,1415,1614]
    ps_n = [7,11,10,7,10,10,7,6,6,4,4,2,2,2,2,2,2]

    vcf_in = vcf.Reader(open('../022117.somatic.filtered.nodbsnp.varscan.vcf','r'))
    primaries = {i:[] for i in ps}

    for variant in vcf_in:
        for s in variant.samples:
            if s.called and sample_map[s.sample] in primaries:
                primaries[sample_map[s.sample]].append(str(variant.CHROM) + '-' + str(variant.POS))

    fout = open('shared_by_primary.csv','w')
    fout.write('count')
    for i in ps:
        fout.write(',' + str(i))
    fout.write('\n')

    for i,j in primaries.items():
        primaries[i] = Counter(Counter(j).values())

    for i in range(1,12):
        fout.write(str(i))
        for p in ps:
            if i in primaries[p]:
                fout.write(',' + str(primaries[p][i]))
            else:
                fout.write(',0')
        fout.write('\n')
    fout.close()

def compute_common_matrix():

    vcf_in = vcf.Reader(open('../022117.somatic.filtered.nodbsnp.varscan.vcf','r'))
    samples = {i:[] for i in vcf_in.samples}

    for variant in vcf_in:
        for s in variant.samples:
            if s.called:
                samples[s.sample].append(str(variant.CHROM) + '-' + str(variant.POS))

    for i,j in samples.items():
        samples[i] = set(j)

    fout = open('common_variants_matrix_count.csv','w')
    fout2 = open('common_variants_matrix_proportion.csv','w')
    fout.write(',' + ','.join(vcf_in.samples) + '\n')
    fout2.write(',' + ','.join(vcf_in.samples) + '\n')
    for i in vcf_in.samples:
        fout.write(i)
        fout2.write(i)
        for j in vcf_in.samples:
            n = len(samples[i].union(samples[j]))

            fout.write(',' + str(len(samples[i].intersection(samples[j]))))
            fout2.write(',' + str(len(samples[i].intersection(samples[j]))/float(n)))
        fout.write('\n')
        fout2.write('\n')
    fout.close()
    fout2.close()


def prepare_for_oncoprint():

    vcf_in = vcf.Reader(open('../022117.somatic.filtered.nodbsnp.varscan.vcf','r'))
    data = []

    for variant in vcf_in:

        for s in variant.samples:
            if s.called:
                temp_effects = []
                effects = []
                for ann in variant.INFO['ANN']:
                    info = ann.split('|')
                    gene = str(info[3])
                    effect = str(info[1])

                    temp = [s.sample,gene,'NULL',effect]
                    temp_effects.append(temp)
                    effects.append(effect)

                if len(temp_effects)==1:
                    data.append(temp_effects[0])
                else:
                    import pdb; pdb.set_trace()

    samples = data.keys()

    fout = open('mutations_oncoprint.txt','w')
    fout.write('Sample\tGene\tAlteration\tType\n')
    for i in samples:
        fout.write(i + ',' + str(len(data[i])))
        for j in ['missense_variant','splice_donor_variant','frameshift_variant','truncating','synonymous']:
            fout.write(',' + str(data[i].count(j)))
        fout.write('\n')

    fout.close()


def prepare_for_plotting():

    sm = SomaticMutations('../somatic.nodbsnp.varscan.filtered.vcf')
    sm.parse_deleterious_types()
    sm.save_deleterious_types('mutation_types.csv')


prepare_for_plotting()
#prepare_for_oncoprint()

#compute_common_matrix()
#compute_shared_within_primary()
#count_by_gene()
