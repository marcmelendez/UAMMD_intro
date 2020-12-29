set term postscript eps enhanced color font "Roman" size 5cm,5cm
set output "flexibleRod.eps"

set size square
unset colorbox
set xlabel "{/:Italic x}"
set ylabel "{/:Italic y}"

set xtics -1,0.5,0.5
set ytics -1,0.5,0.5

set palette defined (0 "grey", 15000 "black")

plot [-0.5:1.0][-1.0:0.5] \
     "../data/cable.dat" u 1:2:0 w l lw 2 lc palette notitle, \
     0 dt 2 lw 2 lc rgb 'black' notitle
