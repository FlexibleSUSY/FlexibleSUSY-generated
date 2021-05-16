set title "MRSSM2 renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='MRSSM2_rgflow.dat'

plot for [i=2:99+1] filename using 1:(column(i)) title columnhead(i)
