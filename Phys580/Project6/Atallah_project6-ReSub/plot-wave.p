# Adrien Atallah 
# Project 6
#
# 'plot-wave.p' plots the evolution of a Gaussian wavefunction as 
# the probability density vs location in the box.  Gaussian
# is initially centered at x = 0, with width = 0.5.  
# Values are taken from 'wave.txt'
  
set terminal postscript
set output 'Atallah_Project6-wave.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "x position"
set ylabel "probability density"

#set xr [-2.0:2.0]
#set yr [-0.5:0.5]


plot  "wave.txt" using 1:2 title 'Evolution of Gaussian Wavefunction' with lines

