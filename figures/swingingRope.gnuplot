set term postscript eps enhanced color font "Roman" size 5cm,5cm
set output "swingingRope.eps"

set size square
unset colorbox
set xlabel "{/:Italic x}"
set ylabel "{/:Italic y}"

set palette defined (0 "grey", 15000 "black")

plot [-1:1][-1.25:0.75] \
     "../data/swingingRope.dat" u 1:2:0 w l lw 2 lc palette notitle, \
     0 dt 2 lw 2 lc rgb 'black' notitle
