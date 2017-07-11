agfusion annotate \
  --gene5prime ENSMUSG00000057841 \
  --gene3prime ENSMUSG00000000441 \
  --junction5prime 115807802 \
  --junction3prime 115626966 \
  --genome GRCm38 \
  --out RPL32-RAF1 \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3

agfusion annotate \
  --gene5prime ENSMUSG00000057841 \
  --gene3prime ENSMUSG00000000441 \
  --junction5prime 115807802 \
  --junction3prime 115626966 \
  --genome GRCm38 \
  --out RPL32-RAF1_nodomain \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3 \
  --no_domain_labels

agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3

agfusion annotate \
  --gene5prime ENSMUSG00000022770 \
  --gene3prime ENSMUSG00000002413 \
  --junction5prime 31684294 \
  --junction3prime 39648486 \
  --genome GRCm38 \
  --out DLG1-BRAF_nodomain \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3 \
  --no_domain_labels

agfusion annotate \
  --gene5prime ENSMUSG00000030849 \
  --gene3prime ENSMUSG00000040265 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  --genome GRCm38 \
  --out FGFR2-DNM3 \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --recolor "TMhelix;black" \
  --rename "Pkinase_Tyr;Kinase" \
  --rename "TMhelix;TM" \
  --scale 2700 \
  --width 20 \
  --height 3 \
  --fontsize 12

agfusion annotate \
  --gene5prime ENSMUSG00000030849 \
  --gene3prime ENSMUSG00000040265 \
  --junction5prime 130167703 \
  --junction3prime 162019992 \
  --genome GRCm38 \
  --out FGFR2-DNM3_nodomain \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --recolor "TMhelix;black" \
  --rename "Pkinase_Tyr;Kinase" \
  --rename "TMhelix;TM" \
  --scale 2700 \
  --width 20 \
  --height 3 \
  --fontsize 12 \
  --no_domain_labels

agfusion annotate \
  --gene5prime ENSMUSG00000030849 \
  --gene3prime ENSMUSG00000055322 \
  --junction5prime 130167703 \
  --junction3prime 74016186 \
  --genome GRCm38 \
  --out FGFR2-TNS1 \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --recolor "TMhelix;black" \
  --rename "Pkinase_Tyr;Kinase" \
  --rename "TMhelix;TM" \
  --width 20 \
  --height 3 \
  --scale 2700 \
  --fontsize 12

agfusion annotate \
  --gene5prime ENSMUSG00000030849 \
  --gene3prime ENSMUSG00000055322 \
  --junction5prime 130167703 \
  --junction3prime 74016186 \
  --genome GRCm38 \
  --out FGFR2-TNS1_nodomain \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --recolor "TMhelix;black" \
  --rename "Pkinase_Tyr;Kinase" \
  --rename "TMhelix;TM" \
  --width 20 \
  --height 3 \
  --scale 2700 \
  --fontsize 12 \
  --no_domain_labels

agfusion annotate \
  --gene5prime ENSMUSG00000042699 \
  --gene3prime ENSMUSG00000000441 \
  --junction5prime 153482742 \
  --junction3prime 115626966 \
  --genome GRCm38 \
  --out DHX9-RAF1 \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3 \
  --fontsize 12

agfusion annotate \
  --gene5prime ENSMUSG00000042699 \
  --gene3prime ENSMUSG00000000441 \
  --junction5prime 153482742 \
  --junction3prime 115626966 \
  --genome GRCm38 \
  --out DHX9-RAF1_nodomain \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3 \
  --fontsize 12 \
  --no_domain_labels

agfusion annotate \
  --gene5prime FGFR2 \
  --gene3prime CCDC6 \
  --junction5prime 123243212 \
  --junction3prime 61612460 \
  --genome GRCh37 \
  --out FGFR2-CCDC6 \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --exclude_domain Pkinase ig \
  --rename "Tmhmm;TM" \
  --width 20 \
  --height 3 \
  --scale 2700 \
  --fontsize 12

agfusion annotate \
  --gene5prime FGFR2 \
  --gene3prime CCDC6 \
  --junction5prime 123243212 \
  --junction3prime 61612460 \
  --genome GRCh37 \
  --out FGFR2-CCDC6_nodomain \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --exclude_domain Pkinase ig \
  --rename "Tmhmm;TM" \
  --width 20 \
  --height 3 \
  --scale 2700 \
  --fontsize 12 \
  --no_domain_labels

agfusion annotate \
  --gene5prime ENSMUSG00000009376 \
  --gene3prime ENSMUSG00000007655 \
  --junction5prime 17513541 \
  --junction3prime 17339113 \
  --genome GRCm38 \
  --out MET-CAV1 \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3

agfusion annotate \
  --gene5prime ENSMUSG00000009376 \
  --gene3prime ENSMUSG00000007655 \
  --junction5prime 17513541 \
  --junction3prime 17339113 \
  --genome GRCm38 \
  --out MET-CAV1_nodomain \
  --middlestar \
  --recolor "Pkinase_Tyr;red" \
  --rename "Pkinase_Tyr;Kinase" \
  --scale 2700 \
  --width 20 \
  --height 3 \
  --no_domain_labels
