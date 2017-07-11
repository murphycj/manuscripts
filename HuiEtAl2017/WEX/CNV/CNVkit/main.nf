params.samtools = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/lib/samtools-1.3.1/samtools'
params.access = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/cnvkit/access-10kb.GRCm38.bed'
params.reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa'

params.brca1_fl_reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/ref/brca1_fl_reference.corrected.cnn'
params.liver_reference = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/ref/livers.cnn'

params.my_targets_sureselect = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/cnvkit/SureSelect/my_targets.bed'
params.my_antitargets_sureselect = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/cnvkit/SureSelect/my_antitargets.bed'

params.my_targets_nimblegen = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/cnvkit/NimbleGen/my_targets.bed'
params.my_antitargets_nimblegen = '/athena/elementolab/scratch/chm2059/from_dat02/chm2059/data/refdata/GRCm38/cnvkit/NimbleGen/my_antitargets.bed'

brca1_fl = Channel.from(
    'HL-WES-01',
    'HL-WES-02',
    'HL-WES-03',
    'HL-WES-04',
    'HL-WES-05',
    'HL-WES-06',
    'HL-WES-07',
    'HL-WES-08',
    'HL-WES-09',
    'HL-WES-10',
    'HL-WES-11'
)


tumors = Channel.from(
    'HL-WES-18',
    'HL-WES-19',
    'HL-WES-20',
    'HL-WES-21',
    'HL-WES-22',
    'HL-WES-23',
    'HL-WES-24',
    'HL-WES-25',
    'HL-WES-26',
    'HL-WES-27',
    'HL-WES-28',
    'HL-WES-29',
    'HL-WES-30',
    'HL-WES-42',
    'HL-WES-43',
    'HL-WES-44',
    'HL-WES-45',
    'HL-WES-46',
    'HL-WES-47',
    'HL-WES-48',
    'HL-WES-49',
    'HL-WES-50',
    'HL-WES-51',
    'HL-WES-52',
    'HL-WES-53',
    'HL-WES-54',
    'HL-WES-55',
    'HL-WES-56',
    'HL-WES-57',
    'HL-WES-58',
    'HL-WES-59',
    'HL-WES-60',
    'HL-WES-61',
    'HL-WES-62',
    'HL-WES-64',
    'HL-WES-65',
    'HL-WES-66',
    'HL-WES-67',
    'HL-WES-68',
    'HL-WES-69',
    'HL-WES-70',
    'HL-WES-71',
    'HL-WES-72',
    'HL-WES-73',
    'HL-WES-74',
    'HL-WES-75',
    'HL-WES-76',
    'HL-WES-77',
    'HL-WES-78',
    'HL-WES-79',
    'HL-WES-80',
    'HL-WES-81'
  )

process cnvkit {

  storeDir "${baseDir}/${tumor}"
  maxForks 30

  input:
    val tumor from tumors

  output:
    set file("*.cnr"), file("*.cns"), file("${tumor}-scatter.png"), file("${tumor}-uncorrected-scatter.png") into cnvkit_get_cnr_out

  """
  export R_HOME=\"\"

  python ~/.local/bin/cnvkit.py fix \
    /athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/coverage/${tumor}.targetcoverage.cnn \
    /athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/coverage/${tumor}.antitargetcoverage.cnn \
    ${params.liver_reference} \
    -o ${tumor}.cnr

  python ~/.local/bin/cnvkit.py segment \
    -m cbs \
    --drop-low-coverage \
    ${tumor}.cnr \
    -o ${tumor}.cns

  python ~/.local/bin/cnvkit.py scatter ${tumor}.cnr -s ${tumor}.cns -o ${tumor}-uncorrected-scatter.png --y-min -4 --y-max 4

  Rscript /athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/correct_variance.R \
    -cnr ${tumor}.cnr \
    -reference ${params.liver_reference} \
    -prefix ${tumor} \
    -ymax 7 \
    -ymin -7 \
    -correction new

  python ~/.local/bin/cnvkit.py segment \
    -m cbs \
    --drop-low-coverage \
    ${tumor}.corrected.cnr \
    -o ${tumor}.corrected.cns

  python /athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/call.py \
    --sample ${tumor} \
    --cns ${tumor}.corrected.cns \
    --cnr ${tumor}.corrected.cnr

  python ~/.local/bin/cnvkit.py scatter ${tumor}.corrected.call.cnr -s ${tumor}.corrected.call.cns -o ${tumor}-scatter.png --y-min -4 --y-max 4
  """
}

process cnvkit_initial_brca1_fl {

  storeDir "${baseDir}/${tumor}"

  input:
    val tumor from brca1_fl

  output:
    set tumor, file('*.cnr'), file('*.cns'), file('*scatter.png'), file('*call.cns'), file('*call.cnr') into cnvkit_get_cnr_brca1_fl_out

  """
  export R_HOME=\"\"

  python ~/.local/bin/cnvkit.py fix \
    /athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/coverage/${tumor}.targetcoverage.cnn \
    /athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/coverage/${tumor}.antitargetcoverage.cnn \
    ${params.brca1_fl_reference} \
    -o ${tumor}.cnr

  python ~/.local/bin/cnvkit.py segment \
    -m cbs \
    ${tumor}.cnr \
    -o ${tumor}.cns

  python /athena/elementolab/scratch/chm2059/from_dat02/chm2059/elementoCantley_2014_9_29/results/WEX/CNV/CNVkit/call.py \
    --sample ${tumor} \
    --cns ${tumor}.cns \
    --cnr ${tumor}.cnr

  python ~/.local/bin/cnvkit.py scatter ${tumor}.call.cnr -s ${tumor}.call.cns -o ${tumor}-scatter.png --y-min -4 --y-max 4
  """
}
