export remaining=~/Desktop/Research/0914_hui/results/WEX/CNV/CNVkit/plots/heatmap
cnvkit.py heatmap \
  $remaining/HL-WES-23.cns \
  $remaining/HL-WES-57.cns \
  $remaining/HL-WES-61.cns \
  $remaining/HL-WES-64.cns \
  $remaining/HL-WES-67.cns \
  $remaining/HL-WES-68.cns \
  $remaining/HL-WES-70.cns \
  $remaining/HL-WES-71.cns \
  $remaining/HL-WES-73.cns \
  $remaining/HL-WES-81.cns \
  -c 6:10000000-30000000 -o chr6.cnv.pdf

cnvkit.py heatmap \
  $remaining/HL-WES-06.cns \
  $remaining/HL-WES-18.cns \
  $remaining/HL-WES-28.cns \
  $remaining/HL-WES-30.cns \
  $remaining/HL-WES-57.cns \
  $remaining/HL-WES-76.cns \
  $remaining/HL-WES-77.cns \
  $remaining/HL-WES-78.cns \
  $remaining/HL-WES-81.cns \
  -c 9:3000000-20000000 -o chr9.cnv.pdf

cnvkit.py heatmap \
  $remaining/HL-WES-49.cns \
  $remaining/HL-WES-51.cns \
  $remaining/HL-WES-72.cns \
  -c 15:55000000-65000000 -o chr15.myc.cnv.pdf

cnvkit.py heatmap \
  $remaining/HL-WES-02.cns \
  $remaining/HL-WES-05.cns \
  $remaining/HL-WES-06.cns \
  $remaining/HL-WES-07.cns \
  $remaining/HL-WES-08.cns \
  $remaining/HL-WES-10.cns \
  $remaining/HL-WES-11.cns \
  $remaining/HL-WES-22.cns \
  $remaining/HL-WES-24.cns \
  $remaining/HL-WES-47.cns \
  $remaining/HL-WES-48.cns \
  $remaining/HL-WES-49.cns \
  $remaining/HL-WES-51.cns \
  $remaining/HL-WES-53.cns \
  $remaining/HL-WES-54.cns \
  $remaining/HL-WES-59.cns \
  $remaining/HL-WES-61.cns \
  $remaining/HL-WES-74.cns \
  -c 14:30000000-120000000 -o chr14.rb1.cnv.pdf

cnvkit.py scatter -s ../HL-WES-09/HL-WES-09.call.cn{s,r} -c 19:30000000-33830000 -g Pten -o HL-WES-09.pten.pdf --title HL-WES-09

cnvkit.py scatter -s ../HL-WES-08/HL-WES-08.call.cn{s,r} -c 11:10000000-20000000 -g Egfr -o HL-WES-08.egfr.pdf --title HL-WES-08

cnvkit.py scatter -s ../HL-WES-04/HL-WES-04.call.cn{s,r} -c 7:128000000-132000000 -g Fgfr2 -o HL-WES-04.fgfr2.pdf --title HL-WES-04

cnvkit.py scatter -s ../HL-WES-22/HL-WES-22.corrected.call.cn{s,r} -c 7:128000000-132000000 -g Fgfr2 -o HL-WES-22.fgfr2.pdf --title HL-WES-22
