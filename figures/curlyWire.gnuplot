set term postscript eps enhanced color font "Roman" size 10cm,5cm
set output "curlyWire.eps"
set xtics -.2,.1,.2
set ytics -.2,.1,.2

set size square

set multiplot layout 1,2

set xlabel "{/:Italic x}"
set ylabel "{/:Italic y}"
plot [-.2:.2][-.2:.2] \
     "< head -n 102 ../data/curlyWire.dat" u 1:2 w l dt 2 lw 3 lc rgb 'black' notitle, \
     "< tail -n 101 ../data/curlyWire.no_torsion.dat" u 1:2 w l lw 3 lc rgb 'grey' notitle, \
     "< tail -n 101 ../data/curlyWire.dat" u 1:2 w l lw 3 lc rgb 'black' notitle
set xlabel "{/:Italic y}"
set ylabel "{/:Italic z}"
plot [-.2:.2][-.2:.2] \
     "< head -n 102 ../data/curlyWire.dat" u 2:3 w l dt 2 lw 3 lc rgb 'black' notitle, \
     "< tail -n 101 ../data/curlyWire.no_torsion.dat" u 2:3 w l lw 3 lc rgb 'grey' notitle, \
     "< tail -n 101 ../data/curlyWire.dat" u 2:3 w l lw 3 lc rgb 'black' notitle

unset multiplot
