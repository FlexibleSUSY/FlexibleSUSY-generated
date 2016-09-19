set terminal x11

set title "MRSSMtower renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='MRSSMtower_rgflow.dat'

plot for [i=2:99+1] filename using 1:(column(i)) title columnhead(i)
