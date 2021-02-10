set term postscript eps enhanced color font "Roman" size 10cm,5cm
set output "rubberBallE.eps"

set multiplot layout 2,1
set xlabel "Time [s]"
set ylabel "Height [m]"
plot [0:7][0:2.5] "../data/rubberBallE.dat" u 1:3 w l lw 3 lc rgb 'black' notitle
set ylabel "Energy [J]"
plot [0:7][0:2.5] "../data/rubberBallE.dat" u 1:4 w l lw 3 lc rgb 'black' notitle
unset multiplot
