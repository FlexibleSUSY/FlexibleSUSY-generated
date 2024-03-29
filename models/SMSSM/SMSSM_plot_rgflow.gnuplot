set title "SMSSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='SMSSM_rgflow.dat'

plot for [i=2:121+1] filename using 1:(column(i)) title columnhead(i)
