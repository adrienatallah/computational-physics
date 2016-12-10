# Adrien Atallah 
# Project 3
#
# 'plot_extra.p' is a Gnuplot script to plot 'fort.10' which contains the temperature 
# vs. the year extrapolated to the year 2050 calculated using the 2nd 
# degree polynomial
  
set terminal postscript
set output 'Atallah_Project3-Extrapolation.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

set title "Extrapolated Temperature using Best Fit Polynomial (2nd order)"
set xlabel "Year"
set ylabel "Temperature"

set xr [1950.0:2052.0]
set yr [50.0:60.0]

plot  "fort.10" using 1:2 title 'Best Fit Polynomial (2nd order)' with lines
