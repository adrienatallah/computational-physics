# Adrien Atallah 
# Project 6
#
# 'plot-01dt.p' plots the width of the Gaussian versus time for
# dt = 0.1
  
set terminal postscript
set output 'Atallah_Project6-01dtvar.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "time"
set ylabel "width"

set label "box: -5 < x < 5" at 0.15,5.25
set label "N=30; dt=0.1" at 0.15,5.0
set label "Initial Gaussian Centered at x = 0" at 0.15,4.75

#set xr [-2.0:2.0]
#set yr [-0.5:0.5]

plot "variance.txt" using 1:2 title 'Evolution of the width of a Gaussian Wavefunction' with lines
