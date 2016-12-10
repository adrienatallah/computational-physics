# Adrien Atallah 
# Project 2
# Gnuplot script to plot error dependence on number of lattice points.
# Data is in 'fort.8' which contains the number of lattice points in the first
# column and the difference between the approximate and exact results in the second.

  
set terminal postscript
set output 'Atallah_Project2-error.ps'                     
set logscale                            
unset label                            
set xtic auto                          
set ytic auto                          

set title "Number of lattice points vs approximation error"
set xlabel "N : Number of lattice points"
set ylabel "Error"
set label "H=10m; W=20m; D=30m; x0=100m; y0=100m" at 8.0,0.000000001




set xr [4.0:100.0]
set yr [0.000000000000000001:0.000001]

plot    "fort.8" using 1:2 title 'Error' with lines
