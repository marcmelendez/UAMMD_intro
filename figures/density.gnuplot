set term postscript eps enhanced color font "Roman" size 10cm,5cm
set output "density.eps"

set xlabel "Height {/:Italic z}"
set ylabel "Density"
set ytics 0,1e-7,5e-7
unset colorbox

set palette defined (0 "grey", 5 "black")
plot [-200:200][0:5e-7] "../data/density.dat" u 1:2:($0/101) w l lw 3 lc palette notitle
