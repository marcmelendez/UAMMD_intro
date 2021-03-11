set term postscript eps enhanced color font "Roman" size 10cm,5cm
set output "Langevin.eps"
set xlabel "Time [{/:Italic m}^{1/2} {/Symbol s} {/Symbol e}^{-1/2}]"

set size square
set key bottom right
set key samplen 5

set multiplot layout 1,2
set arrow 1 from 10,0 to 10,2.5 lw 2 dt 3 lc rgb 'gray70' nohead
set ylabel "Energy per particle [{/Symbol e}]"
plot [0:100][0:2.5] \
     "../data/LJmacro.NVE.dat" u 1:($2/1e5) w l lw 4 lc rgb 'gray60' title "NVE", \
     "../data/LJmacro.NVT.friction_0.1.dat" u ($1+10):($2/1e5) w l lw 3 dt 1 lc rgb 'black' title "NVT", \
     "../data/LJmacro.NVT.friction_0.05.dat" u ($1+10):($2/1e5) w l lw 3 dt 4 lc rgb 'black' notitle, \
     "../data/LJmacro.NVT.friction_0.025.dat" u ($1+10):($2/1e5) w l lw 3 dt 5 lc rgb 'black' notitle
unset arrow 1

set arrow 2 from 10,0 to 10,1.75 lw 2 dt 3 lc rgb 'gray70' nohead
set ylabel "Thermal energy [{/Symbol e}]"
plot [0:100][0:1.75] \
     "" w l lw 3 dt 2 lc rgb 'gray60' title "Target", \
     "../data/LJmacro.NVE.dat" u 1:6 w l lw 4 lc rgb 'gray60' notitle, \
     "../data/LJmacro.NVT.friction_0.1.dat" u ($1+10):6 w l lw 3 dt 1 lc rgb 'black' title "{/Symbol g} = 0.1", \
     "../data/LJmacro.NVT.friction_0.05.dat" u ($1+10):6 w l lw 3 dt 4 lc rgb 'black' title "{/Symbol g} = 0.05", \
     "../data/LJmacro.NVT.friction_0.025.dat" u ($1+10):6 w l lw 3 dt 5 lc rgb 'black' title "{/Symbol g} = 0.025", \
     1.5 w l lw 3 dt 2 lc rgb 'gray60' notitle

unset multiplot
