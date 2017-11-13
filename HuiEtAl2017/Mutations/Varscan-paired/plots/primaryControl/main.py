from collections import Counter
import vcf
import pandas
from seqpy.parsers import *

def prepare_for_plotting():

    aa = SnpEffAnnotation(classes=[ 'missense','synonymous','indel','frameshift','nonsense','splicesite'])

    sm = SummarizeVCFAnnotation(filename='../../somatic.nodbsnp.varscan.filtered.primaryControl.vcf',annotation=aa)
    sm.summarize()
    sm.save_summary('mutation_types.csv')


prepare_for_plotting()
