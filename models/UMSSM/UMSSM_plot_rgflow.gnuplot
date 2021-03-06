set title "UMSSM renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='UMSSM_rgflow.dat'

plot for [i=2:142+1] filename using 1:(column(i)) title columnhead(i)
