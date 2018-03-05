set title "MSSMNoFVHimalaya renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='MSSMNoFVHimalaya_rgflow.dat'

plot for [i=2:111+1] filename using 1:(column(i)) title columnhead(i)
