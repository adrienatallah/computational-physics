# Adrien Atallah 
# Gnuplot script to plot project 1 data
# plotting "fort.8" which includes dx vs the difference between the 3-pt formula and exact 2nd derivative
# also plotting "fort.9" which includes dx vs the difference between the 5-pt formula and exact 2nd derivative

#set   autoscale   
set output 'Atallah_Project1.eps'                     
set logscale                            
unset label                            # remove any previous labels
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically

set title "dx vs. Difference between 3 and 5 pt formula approx. for 2nd derivative of xsinx and exact"
set xlabel "dx"
set ylabel "approx - exact"



set xr [0.0001:0.1]
set yr [0.001:5]

plot    "fort.8" using 1:2 title '3pt' with lines,\
	"fort.9" using 1:2 title '5pt' with lines
