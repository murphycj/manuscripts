python /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/chm2059/lib/python2.7.8/bin/cnvkit.py \
	access ../../Mus_musculus.GRCm38.dna_sm.primary_assembly.fa \
	-s 10000 \
	-o access-10kb.GRCm38.bed
python /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/chm2059/lib/python2.7.8/bin/cnvkit.py target ~/chm2059/data/Mouse_Exome_Design/110624_MM10_exome_L2R_D02_EZ_HX1-v4.bed --split --annotate ../../082016.refFlat_ensembl.txt -o my_targets.bed
python /zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/chm2059/lib/python2.7.8/bin/cnvkit.py antitarget ~/chm2059/data/Mouse_Exome_Design/110624_MM10_exome_L2R_D02_EZ_HX1-v4.bed --access access-10kb.GRCm38.bed -o my_antitargets.bed
