# Adrien Atallah 
# Project 5
#
# 'plotb2.p' plots an elliptical orbit with values as if the earth where 
# orbiting the sun.  Values from 'locations.txt' when using 'proj5basic2.input'
# with 'proj5.script'.
  
set terminal postscript
set output 'Atallah_Project5-OrbitsB2.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "x position (meters)"
set ylabel "y position (meters)"
set label "x1=1.0E11m; y1=1.0E11m" at -1.49E11,1.4E11
set label "vx1=-20.0E3m/s; vy1=20.0E3m/s" at -1.49E11,1.3E11
set label "dt=1.5E5s" at -1.49E11,1.2E11


#set xr [-2.0:2.0]
#set yr [-0.5:0.5]


plot  "locations.txt" using 2:3 title 'Orbit of planet 1 with mass m1' with lines
