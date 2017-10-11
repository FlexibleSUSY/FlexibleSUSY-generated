set title "HTHDMIIMSSMBC renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='HTHDMIIMSSMBC_rgflow.dat'

plot for [i=2:43+1] filename using 1:(column(i)) title columnhead(i)
