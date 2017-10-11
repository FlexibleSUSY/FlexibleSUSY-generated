set title "SplitMSSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='SplitMSSM_rgflow.dat'

plot for [i=2:41+1] filename using 1:(column(i)) title columnhead(i)
