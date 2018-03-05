set title "THDMII renormalization group flow"
set xlabel "renormalization scale / GeV"
set logscale x

if (!exists("filename")) filename='THDMII_rgflow.dat'

plot for [i=2:42+1] filename using 1:(column(i)) title columnhead(i)
