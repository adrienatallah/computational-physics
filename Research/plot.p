# Adrien Atallah 
#
# 'plot.p' plots ground state and first excited state wave function for 
# a harmonic oscillator.
  
set terminal postscript
set output 'Harmonic Oscillator.ps'                     
unset logscale                            
unset label 
set autoscale                           
set xtic auto                          
set ytic auto                          

unset title 
set xlabel "x"
set ylabel "Wave Function - Psi(x)"

#set xr [-2.0:2.0]
#set yr [-0.5:0.5]


plot  "GroundStateHO.txt" using 1:2 title 'Harmonic Oscillator - Ground State' with lines
plot  "1stExStateHO.txt" using 1:2 title 'Harmonic Oscillator - First Excited State' with lines
plot  "2ndExStateHO.txt" using 1:2 title 'Harmonic Oscillator - Second Excited State' with lines
plot  "3rdExStateHO.txt" using 1:2 title 'Harmonic Oscillator - Third Excited State' with lines
plot  "4thExStateHO.txt" using 1:2 title 'Harmonic Oscillator - Fourth Excited State' with lines
plot  "5thExStateHO.txt" using 1:2 title 'Harmonic Oscillator - Fifth Excited State' with lines
