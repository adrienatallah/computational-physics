# Adrien Atallah 
# Project 5
#
# 'plota1.p' plots the positions of 2 bodies orbiting a central primary mass. 
# Values from 'locations.txt' when using 'proj5advanced1.input'
# with 'proj5.script'.
  
set terminal postscript
set output 'Atallah_Project5-OrbitsA1.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "x position (meters)"
set ylabel "y position (meters)"
set label "x1=1.0E11m; y1=1.0E11m" at -1.9E11,1.9E11
set label "x2=0.0m; y2=1.5E11m" at -1.9E11,1.75E11
set label "vx1=-20.3E3m/s; vy1=20.3E3m/s" at -1.9E11,1.6E11
set label "vx2=-3.0E4E3m/s; vy2=0.0m/s" at -1.9E11,1.45E11
set label "dt=1.5E5s" at -1.9E11,1.3E11
 

plot    "locations.txt" using 2:3 title '1st body (m1) orbiting M' with lines,\
	"locations.txt" using 4:5 title '2nd body (m2) orbiting M' with lines
