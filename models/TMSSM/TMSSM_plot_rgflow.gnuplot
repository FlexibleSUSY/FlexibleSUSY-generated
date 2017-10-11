set title "TMSSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='TMSSM_rgflow.dat'

plot for [i=2:117+1] filename using 1:(column(i)) title columnhead(i)
