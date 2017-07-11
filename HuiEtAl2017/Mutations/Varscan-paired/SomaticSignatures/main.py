
from collections import Counter
import vcf
import pandas

samples = pandas.read_excel("/Users/charlesmurphy/Desktop/Research/0914_hui/data/samples.xlsx","samples")

brca1 = dict(
    zip(samples['sample_ID'],samples['brca1'])
)

mouse = dict(
    zip(samples['sample_ID'],samples['mouse'])
)


vcf_in = vcf.Reader(open('../somatic.nodbsnp.varscan.filtered.vcf','r'))
data = {i:[] for i in vcf_in.samples}

fout = open('variant_table.csv','w')
fout.write('chrom,pos,ref,alt,sample,brca1,mouse\n')

for variant in vcf_in:

    try:
        #fitler out the indels
        if len(variant.ALT[0])>1 or len(variant.REF)>1 or len(variant.ALT)>1:
            continue
    except:
        import pdb; pdb.set_trace()

    chrom = 'chr' + str(variant.CHROM)

    for s in variant.samples:
        if s.called:
            if brca1[s.sample]=='flox/flox':
                brca='DEL'
            else:
                brca="WT"
            fout.write(
                str(chrom) + ',' +
                str(variant.POS) + ',' +
                str(variant.REF) + ',' +
                str(variant.ALT[0]) + ',' +
                str(s.sample) + ',' +
                str(brca) + ',' +
                'P.' + str(mouse[s.sample]) + '\n'
                )

fout.close()
