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

set output "output/plot5.eps"

set multiplot layout 4,1 rowsfirst

set xlabel '{/:Bold x}' font ",18"
set ylabel '{/:Bold Density}' font ",18"
#set yr[0:1.2]
set key top left
plot "output/txt/test5/field1MUSCL64noAMR.txt" using 1:2 title "N=64, no AMR" with points pointtype 7 pointsize 0.6 linecolor rgb "#3355FF"  lw 6 lt 5,\
"output/txt/test5/field1MUSCL64AMR.txt" using 1:2 title "N=64, AMR" with points pointtype 9 pointsize 0.6 linecolor rgb "#FF0000"  lw 6 lt 5,\
"output/txt/test5/field1MUSCLexact.txt" using 1:2 title "Exact Riemann" with lines linecolor rgb "#000000"  lw 4 lt 3

set xlabel '{/:Bold x}' font ",18"
set ylabel '{/:Bold Velocity}' font ",18"
#set yr[0:1.2]
plot "output/txt/test5/field1MUSCL64noAMR.txt" using 1:3 title "" with points pointtype 7 pointsize 0.6 linecolor rgb "#3355FF"  lw 6 lt 5,\
"output/txt/test5/field1MUSCL64AMR.txt" using 1:3 title "" with points pointtype 9 pointsize 0.6 linecolor rgb "#FF0000"  lw 6 lt 5,\
"output/txt/test5/field1MUSCLexact.txt" using 1:3 title "" with lines linecolor rgb "#000000"  lw 4 lt 3

set xlabel '{/:Bold x}' font ",18"
set ylabel '{/:Bold Pressure}' font ",18"
#set yr[0:1.2]
plot "output/txt/test5/field1MUSCL64noAMR.txt" using 1:4 title "" with points pointtype 7 pointsize 0.6 linecolor rgb "#3355FF"  lw 6 lt 5,\
"output/txt/test5/field1MUSCL64AMR.txt" using 1:4 title "" with points pointtype 9 pointsize 0.6 linecolor rgb "#FF0000"  lw 6 lt 5,\
"output/txt/test5/field1MUSCLexact.txt" using 1:4 title "" with lines linecolor rgb "#000000"  lw 4 lt 3


set xlabel '{/:Bold x}' font ",18"
set ylabel '{/:Bold Specific internal energy}' font ",18"
#set yr[*:*]
plot "output/txt/test5/field1MUSCL64noAMR.txt" using 1:5 title "" with points pointtype 7 pointsize 0.6 linecolor rgb "#3355FF"  lw 6 lt 5,\
"output/txt/test5/field1MUSCL64AMR.txt" using 1:5 title "" with points pointtype 9 pointsize 0.6 linecolor rgb "#FF0000"  lw 6 lt 5,\
"output/txt/test5/field1MUSCLexact.txt" using 1:5 title "" with lines linecolor rgb "#000000"  lw 4 lt 3



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


