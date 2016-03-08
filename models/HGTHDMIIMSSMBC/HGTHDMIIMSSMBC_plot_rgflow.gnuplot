set terminal x11

set title "HGTHDMIIMSSMBC renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='HGTHDMIIMSSMBC_rgflow.dat'

plot for [i=2:50+1] filename using 1:(column(i)) title columnhead(i)
