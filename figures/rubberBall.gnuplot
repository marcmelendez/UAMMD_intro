set term postscript eps enhanced color font "Roman" size 5cm,5cm
set output "rubberBall.eps"

set size square;
set xlabel "{/:Italic x}"
set ylabel "{/:Italic y}"

plot [0:2.25][0:2.25] "../data/rubberBall.dat" u 2:3 w p pt 7 ps 1.2 lc rgb 'black' notitle
