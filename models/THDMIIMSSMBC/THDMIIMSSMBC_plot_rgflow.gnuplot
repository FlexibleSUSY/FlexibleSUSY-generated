set title "THDMIIMSSMBC renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='THDMIIMSSMBC_rgflow.dat'

plot for [i=2:42+1] filename using 1:(column(i)) title columnhead(i)
