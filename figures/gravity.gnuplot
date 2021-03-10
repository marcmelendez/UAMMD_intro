set term postscript eps enhanced color font "Roman" size 10cm,4cm
set output "gravity.eps"

set size ratio .667
set xlabel "{/:Italic x} [a.u.]"
set ylabel "{/:Italic y} [a.u.]"
unset colorbox

set palette defined (0 "grey", 5 "black")
set multiplot
plot [-20:40][-20:20] \
     "../data/gravity.dat" u 1:2:($0/21) w p pt 7 ps .1 lc palette notitle, \
     "< tail -n 21 ../data/gravity.dat" u 1:2 w p pt 7 ps 1 lc rgb 'black' notitle
set label 1 "Pluto" at 27,-15
set label 2 "Uranus" at -15,7
set label 3 "Saturn" at 5,-10
set label 4 "Jupiter" at 5,5
plot [-20:40][-20:20] 1/0 notitle
unset multiplot
