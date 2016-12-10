# Adrien Atallah 
# Project 5
#
# 'plot.p' plots ground state and first excited state wave function for 
# the particle in a box, the harmonic oscillator and the Woods-Saxon potential
  
set terminal postscript
set output 'Atallah_Project5-Orbits.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "x position"
set ylabel "y position"

#set xr [-2.0:2.0]
#set yr [-0.5:0.5]


plot  "locations.txt" using 2:3 title 'Orbit of planet 1 with mass m1' with lines
plot  "locations.txt" using 4:5 title 'Orbit of planet 2 with mass m2' with lines
