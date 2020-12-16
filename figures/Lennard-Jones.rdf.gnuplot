set term postscript eps enhanced color font "Roman" size 10cm,5cm
set output "Lennard-Jones.rdf.eps"
set xlabel "{/:Italic r} / {/Symbol s}"
set ylabel "Radial distribution function, {/:Italic g}({/:Italic r})" offset 2,0

set size square

set multiplot layout 1,2

set title "Ideal gas"
plot [0:5][0:2.75] \
     "ideal_gas.rdf.dat" u 1:2:3 w errorbars lc rgb 'black' notitle, \
     1 w l lw 3 lc 'black' dt 2 notitle
set title "Lennard-Jones fluid"
plot [0:5][0:2.75] \
     "Lennard-Jones.rdf.dat" u 1:2:3 w errorlines lw 2 lc rgb 'black' notitle, \
     1 w l lw 3 lc 'black' dt 2 notitle

unset multiplot
