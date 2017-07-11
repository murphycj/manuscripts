import pyensembl
import pandas
import numpy as np

fusions = pandas.read_csv('fusions-noNormal.noReadthrough.n2.csv')

homology = pandas.read_csv("../Drivers/mart_export.txt")
homology = homology[~pandas.isnull(homology['Mouse gene stable ID'])]
homologs = dict(zip(homology['Mouse gene stable ID'],homology['Gene stable ID']))

cosmic_census = pandas.read_csv('/Users/charlesmurphy/Desktop/Research/data/COSMIC/Census_allThu_Jun_8_23_48_12_2017.csv')
cosmic_census_genes = cosmic_census['Gene Symbol'].tolist()

human_db = pyensembl.EnsemblRelease(84,'human')

in_cosmic_5prime = []
in_cosmic_3prime = []

for g in fusions['Gene_1_id.5end_fusion_partner.']:
    if g in homologs:
        gene = human_db.gene_by_id(homologs[g])
        if gene.name in cosmic_census_genes:
            in_cosmic_5prime.append(True)
        else:
            in_cosmic_5prime.append(False)
    else:
        in_cosmic_5prime.append(False)

for g in fusions['Gene_2_id.3end_fusion_partner.']:
    if g in homologs:
        gene = human_db.gene_by_id(homologs[g])
        if gene.name in cosmic_census_genes:
            in_cosmic_3prime.append(True)
        else:
            in_cosmic_3prime.append(False)
    else:
        in_cosmic_3prime.append(False)

in_cosmic_5prime = np.array(in_cosmic_5prime)
in_cosmic_3prime = np.array(in_cosmic_3prime)

fusions = fusions[in_cosmic_3prime | in_cosmic_5prime]
fusions.to_csv('fusions-inCosmicCensus.csv')
