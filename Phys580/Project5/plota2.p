# Adrien Atallah 
# Project 5
#
# 'plota2.p' plots the position of the moon orbiting the Earth,
# treating the Earth as the primary and the moon as m1. 
# the Values are from 'locations.txt' 
# when using 'proj5advanced2.input' with 'proj5.script'.
  
set terminal postscript
set output 'Atallah_Project5-OrbitsA2.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "x position (meters)"
set ylabel "y position (meters)"
set label "x1=0.0m; y1=3.85E8m" at -3.9E8,3.8E8
set label "vx1=-1.023E3m/s; vy1=0.0m/s" at -3.9E8,3.5E8
set label "dt=0.5E5s" at -3.9E8,3.2E8
   

plot    "locations.txt" using 2:3 title 'Moon orbiting the Earth' with lines
