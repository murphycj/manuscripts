### External data preprocessing
```
Data/CaptureKits
```
Data pre-processing involving the two whole exome capture kits.

```
Data/GDAC
```
Data pre-processing involving the TCGA BRCA data downloaded from GDAC.

```
Data/GRCm38
```
Preparation of CNVkit reference data and liftOver of the RNA-editing sites

### Whole-exome sequencing
```
WEX/
```
The whole-exome-seq pipeline for QC, mapping, and initial mutation calling

```
WEX/CNV/CNVkit
```
The copy number analysis pipeline.

### RNAseq

```
Fusions/
```
Post-processing and plotting for the gene fusions.

```
RNAseq/
```
The RNA-seq pipeline for QC, mapping, quantifying, fusions, and initial mutation calling.

```
RNAseq/Cufflinks/
```
Post-processing steps for gene expression estimates (e.g. convert mouse gene symbols to human gene symbols).

```
RNAseq/SubtypeClassification/primaryControl/AIM_signature
```
Scripts to classify tumors with AIMs.

```
RNAseq/SubtypeClassification/primaryControl/PAM50
```
Scripts to classify tumors with PAM50.

```
RNAseq/SubtypeClassification/primaryControl/Markers
```
Scripts to classify tumors according to ER, PR, and HER2 status.

### Mutations
```
Mutations/CountTriNucleotideCoverage
```

Coverage calcuations used for somatic signatures.

```
Mutations/Varscan-paired
```

Mutation calling.

```
Mutations/Varscan-paired/plots/primaryControl
```

Plots related to mutation calling.

```
Mutations/Varscan-paired/excluded_genes.txt
```

List of genes where if any mutation falls within them is removed.

```
Mutations/Varscan-paired/SomaticSignatures/primaryControl
```

Scripts to plot mutational signatures of mouse tumors.