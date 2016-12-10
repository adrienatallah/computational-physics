# Adrien Atallah 
# Project 4
# Adrien Atallah
# 'plotWS8R.p' plots the ground state, 1st, 2nd and 3rd excited states for 
# the Woods-Saxon Potential for R = 8fm
  
set terminal postscript
set output 'Atallah_Project4-WS_R8.ps'
         
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "x"
set ylabel "Wave Function - Psi(x)"

#set xr [-2.0:2.0]
#set yr [-0.5:0.5]


plot  "GroundStateWS.txt" using 1:2 title 'Woods-Saxon Potential - Ground State (R=8fm)' with lines
plot  "1stExStateWS.txt" using 1:2 title 'Woods-Saxon Potential - First Excited State (R=8fm)' with lines
plot  "2ndExStateWS.txt" using 1:2 title 'Woods-Saxon Potential - Second Excited State (R=8fm)' with lines
plot  "3rdExStateWS.txt" using 1:2 title 'Woods-Saxon Potential - Third Excited State (R=8fm)' with lines
