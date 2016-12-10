# Adrien Atallah 
# Project 3
#
# 'plot.p' is a Gnuplot script to plot 'dataset4.dat' with 'fort.7' and 'fort.8'.
# 'dataset4.dat' contains 15 temperature measurements organized
# into 3 columns; the first being the year, the second being the temperature measured
# that year and the third column contains the uncertainty on that measurement.
# 'fort.7' is the temperature vs. the year calculated by the 2nd degree polynomial and 
# 'fort.8' is the same data calculated by a 3rd degree polynomial

  
set terminal postscript
set output 'Atallah_Project3-BestFit.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

set title "Best fit polynomials and data with error bars"
set xlabel "Year"
set ylabel "Temperature"
#set label "T = c1 + c2*y + c3*y^2" at 1932.0, 58.0

set xr [1900.0:2000.0]
set yr [50.0:60.0]

plot    "fort.7" using 1:2 title 'Best Fit Polynomial (2nd order)' with lines,\
	"fort.8" using 1:2 title 'Best Fit Polynomial (3rd order)' with lines,\
	"dataset4.dat" using 1:2:3 title 'Dataset 4' with yerrorbars
