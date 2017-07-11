python /Users/charlesmurphy/Desktop/Research/lib/seqpy/bin/bed_to_kmer_frequency.py \
  --fasta ~/Desktop/Research/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa
  --kmer 3 \
  --bed 110624_MM10_exome_L2R_D02_EZ_HX1-v3.bed \
  --out GRCm38.110624bed.3merCount.csv

python /Users/charlesmurphy/Desktop/Research/lib/seqpy/bin/bed_to_kmer_frequency.py \
  --fasta ~/Desktop/Research/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa
  --kmer 3 \
  --bed S0276129_Regions-mm10-v2.bed \
  --out GRCm38.S0276129.bed.3merCount.csv


cat 110624_MM10_exome_L2R_D02_EZ_HX1-v3.bed S0276129_Regions-mm10-v2.bed | ~/Desktop/tools/bedtools2/bin/bedtools sort -i stdin -faidx names.txt | ~/Desktop/tools/bedtools2/bin/bedtools merge -i stdin > 110624_S0276129.bed

python /Users/charlesmurphy/Desktop/Research/lib/seqpy/bin/bed_to_kmer_frequency.py \
  --fasta ~/Desktop/Research/data/refdata/GRCm38/Mus_musculus.GRCm38.dna.primary_assembly.fa \
  --kmer 3 \
  --bed 110624_S0276129.bed \
  --out GRCm38.110624_S0276129.3merCount.csv
