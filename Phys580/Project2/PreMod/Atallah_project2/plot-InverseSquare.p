# Adrien Atallah 
# Project 2
# Gnuplot script to plot flux from the reactor, flux from a point source
# and the flux using Monte Carlo Integration with error bars
# while the distance between source and detector is varied.
# Data is in 'fort.9' which contains the value of x0 in the first column, 
# the flux from the reactor in the second column, the flux from a point source 
# in the third, the flux using Monte Carlo Integration in the fourth
# and the error of the Monte Carlo approximation in the fifth.

  
set terminal postscript
set output 'Atallah_Project2-InverseSquare.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

set title "Inverse square law: relation between reactor, point source and monte carlo approx. vs x0"
set xlabel "x0: distance between source and detector"
set ylabel "Flux"
set label "H=10m; W=20m; D=30m; y0=10m; N=40" at 30.0,4.0

set xr [1.0:100.0]
set yr [0.001:7.0]

plot    "fort.9" using 1:2 title 'Reactor Flux' with lines,\
	"fort.9" using 1:3 title 'Point Source Flux' with lines,\
	"fort.9" using 1:4:5 title 'Monte Carlo Integration' with yerrorbars
