#!/bin/bash
#This is the plot template of absoption spectrum for different sigma.
Sigma=$(seq 0.001 0.002 0.069)
for sigma in $Sigma
do
	cat >spectrum_Varysigma.gnu <<EOF
set encoding iso_8859_1
# set terminal  postscript enhanced color font "TimesNewRoman, 11" size 5, 4
set terminal  pngcairo truecolor enhanced lw 5.0 font "TimesNewRoman, 44" size 1920, 1680
#set palette rgbformulae 22, 13, -31
# set palette rgbformulae 7,5,15
set key autotitle columnhead
set output "absorption_sigma$sigma.png"
set border
#unset colorbox
set title "absorption\\\_sigma=$sigma" offset 0, -0.8 font "TimesNewRoman, 54"
#set style data linespoints
unset ztics
set key
# set key outside top vertical center
# set pointsize 0.3
set view 0,0
set xtics font "TimesNewRoman, 44"
set xtics offset 0, 0.3
set ytics font "TimesNewRoman, 40"
set xtics 100
set mxtics 2
#set ytics -4, 2, 4
set ylabel font "TimesNewRoman, 48"
set ylabel offset 1.0, 0
set xrange [1100:1700]
set xlabel "wavelength(nm)"
set ylabel "absorption"
#set yrange [-4:4]
#set xtics ("M" 0.00000, "G" 0.888, "K" 1.913, "M" 2.426)
plot './data/spectrum_continuous_Er3_sigma$sigma.dat' using 2:(log10(\$3)) title 'Er3 doped' with lines lc "black",'./data/spectrum_continuous_Bi3-Er3_sigma$sigma.dat' using 2:(log10(\$3)) title 'Bi3/Er3 codoped' with lines lc "red"

EOF
         gnuplot spectrum_Varysigma.gnu
	 mv absorption_sigma$sigma.png figure/

done

