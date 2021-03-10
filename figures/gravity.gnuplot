set term postscript eps enhanced color font "Roman" size 6cm,4cm
set output "gravity.eps"
set palette defined (0 "grey", 5 "black")

set size ratio .667
set xlabel "{/:Italic x} [a.u.]"
set ylabel "{/:Italic y} [a.u.]"
unset colorbox

set label 1 "Pluto" at 27,-15 front
set label 2 "Uranus" at -15,7 front
set label 3 "Saturn" at 5,-10 front
set label 4 "Jupiter" at 5,5 front
plot [-20:40][-20:20] \
     "../data/gravity.dat" u 1:2:($0/21) w p pt 7 ps .1 lc palette notitle, \
     "< tail -n 21 ../data/gravity.dat" u 1:2 w p pt 7 ps 1 lc rgb 'black' notitle
