########################################################################
#                                                                      #
#                             GNUPLOT                                  #
#                                                                      #
########################################################################


# 1D PLOT 

set terminal postscript portrait enhanced color
set fit quiet
set size 3.0/3.0, 4.0/4.0
set size ratio 0.4/1

set output "output/plot1.eps"

set multiplot layout 4,1 rowsfirst

set xlabel '{/:Bold x}' font ",18"
set ylabel '{/:Bold Pressure}' font ",18"
#set yr[0:1.2]
set key top left
plot for [i=1:10] 'output/txt/test8/1MUSCL128time'.i.'.txt' using 1:4 title "" with lines linecolor rgb "#000000"  lw 3 lt 3

set xlabel '{/:Bold x}' font ",18"
set ylabel '{/:Bold Temperature}' font ",18"
#set yr[0:1.2]
plot for [i=1:10] 'output/txt/test8/1MUSCL128time'.i.'.txt' using 1:7 title "" with lines linecolor rgb "#000000"  lw 3 lt 3


set xlabel '{/:Bold x}' font ",18"
set ylabel '{/:Bold Unburnt mass fraction}' font ",18"
#set yr[0:1.2]
plot for [i=1:10] 'output/txt/test8/1MUSCL128time'.i.'.txt' using 1:6 title "" with lines linecolor rgb "#000000"  lw 3 lt 3


set xr[0:1]
set yr[*:*]
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid ls 12
set style line 11 lc rgb '#808080' lt 1
set border 3 ls 11
set key box
set key tmargin
set key font ",12
set key samplen 4  spacing 1.25 
#set output "output/toro1.eps"
#replot 
#set term x11


