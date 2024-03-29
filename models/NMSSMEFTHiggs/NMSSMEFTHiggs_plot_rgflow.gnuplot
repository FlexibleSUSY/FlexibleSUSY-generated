set title "NMSSMEFTHiggs renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='NMSSMEFTHiggs_rgflow.dat'

plot for [i=2:115+1] filename using 1:(column(i)) title columnhead(i)
