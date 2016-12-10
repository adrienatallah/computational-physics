# Adrien Atallah 
# Project 4
#
# 'plot.p' plots ground state and first excited state wave function for 
# the particle in a box, the harmonic oscillator and the Woods-Saxon potential
  
set terminal postscript
set output 'Atallah_Project4-WaveFunction.ps'                     
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


plot  "GroundState.txt" using 1:2 title 'Particle in a Box - Ground State' with lines
plot  "1stExState.txt" using 1:2 title 'Particle in a Box - First Excited State' with lines
plot  "GroundStateHO.txt" using 1:2 title 'Harmonic Oscillator - Ground State' with lines
plot  "1stExStateHO.txt" using 1:2 title 'Harmonic Oscillator - First Excited State' with lines
plot  "GroundStateWS.txt" using 1:2 title 'Woods-Saxon Potential - Ground State' with lines
plot  "1stExStateWS.txt" using 1:2 title 'Woods-Saxon Potential - First Excited State' with lines
