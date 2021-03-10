set term postscript eps enhanced color font "Roman" size 5cm,4cm
set output "Jupiter.eps"

set size square
set xlabel "{/:Italic x} [au]"
set ylabel "{/:Italic y} [au]"
set ytics -3.4,0.1,-3.1
set xtics 3.6,0.1,4
unset colorbox

set palette defined (0 "grey", 5 "black")
set label 1 "Jupiter" at 3.85,-3.12
plot [3.6:4.0][-3.45:-3.05] \
     "< tail -n 10000 ../data/Jupiter.dat" u 1:2:0 w p pt 7 ps 0.1 lc palette notitle, \
     "< tail -n 21 ../data/Jupiter.dat | head -n 10" u 1:2 w p pt 7 ps 1 lc rgb 'black' notitle
