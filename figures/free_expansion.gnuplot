set term postscript eps enhanced color font "Roman" size 10cm,5cm
set output "free_expansion.eps"
set xlabel "{/:Italic x}"
set ylabel "{/:Italic z}"
set xtics -80,40,80
set ytics -80,40,80

set size square

set multiplot layout 1,2
plot [-96:96][-96:96] "< tail -n 100000 ../data/free_expansion.diff.dat" u 1:3 w d lc rgb 'black' notitle, \
                      "< tail -n 100000 ../data/free_expansion.diff.dat" u ($1+128):3 w d lc rgb '0x888888' notitle, \
                      "< tail -n 100000 ../data/free_expansion.diff.dat" u ($1-128):3 w d lc rgb '0x888888' notitle
plot [-96:96][-96:96] "< tail -n 100000 ../data/free_expansion.box.dat" u 1:3 w d lc rgb 'black' notitle, \
                      "< tail -n 100000 ../data/free_expansion.box.dat" u ($1+128):3 w d lc rgb '0x888888' notitle, \
                      "< tail -n 100000 ../data/free_expansion.box.dat" u ($1-128):3 w d lc rgb '0x888888' notitle
unset multiplot
