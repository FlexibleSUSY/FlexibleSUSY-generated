set terminal x11

set title "MSSMRHN renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='MSSMRHN_rgflow.dat'

plot for [i=2:156+1] filename using 1:(column(i)) title columnhead(i)
