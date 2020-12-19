set term postscript eps enhanced color font "Roman" size 5cm,5cm
set output "vibratingString.eps"

set size square
set ytics ("n = 1" 0.2, "n = 2" 0.4, "n = 3" 0.6, "n = 4" 0.8)

plot [0:1][0:1] "../data/vibratingString.n1.dat" u 1:($2 + .2) w l lw 2 lc rgb 'black' notitle, \
                "../data/vibratingString.n2.dat" u 1:($2 + .4) w l lw 2 lc rgb 'black' notitle, \
                "../data/vibratingString.n3.dat" u 1:($2 + .6) w l lw 2 lc rgb 'black' notitle, \
                "../data/vibratingString.n4.dat" u 1:($2 + .8) w l lw 2 lc rgb 'black' notitle
