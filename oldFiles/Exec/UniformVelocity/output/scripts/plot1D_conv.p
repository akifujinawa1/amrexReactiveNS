########################################################################
#                                                                      #
#                             GNUPLOT                                  #
#                                                                      #
########################################################################


# 1D PLOT 

set terminal postscript landscape enhanced color
set fit quiet
set size 3.0/3.0, 4.0/4.0
set size ratio 0.3/0.7

set output "output/plot2.eps"

set multiplot layout 3,2 rowsfirst

set logscale xy
set key font ",10"
set format xy "10^{%+1T}"

set xlabel '{/:Bold n_{Cells}}' font ",18"
set ylabel '{/:Bold E(t_{End})}' font ",18"
#set yr[0:1.2]
set key bottom left
plot "output/txt/conv/toro1error.txt" using 2:3 title "no AMR (1)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#FF0000" lw 2 lt 2,\
"output/txt/conv/toro1error.txt" using 2:8 title "AMR (1)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#000000" lw 2 lt 2,\
"output/txt/conv/toro1error.txt" using 2:6 title "{/Symbol a}=1" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#0D996A"  lw 2 lt 2,\
"output/txt/conv/toro1error.txt" using 2:7 title "{/Symbol a}=2" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#3355FF"  lw 2 lt 2

set xlabel '{/:Bold n_{Cells}}' font ",18"
set ylabel '{/:Bold E(t_{End})}' font ",18"
#set yr[0:1.2]
set key bottom left
plot "output/txt/conv/toro2error.txt" using 2:3 title "no AMR (2)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#FF0000" lw 2 lt 2,\
"output/txt/conv/toro2error.txt" using 2:8 title "AMR (2)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#000000" lw 2 lt 2,\
"output/txt/conv/toro2error.txt" using 2:6 title "{/Symbol a}=1" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#0D996A"  lw 2 lt 2,\
"output/txt/conv/toro2error.txt" using 2:7 title "{/Symbol a}=2" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#3355FF"  lw 2 lt 2

set xlabel '{/:Bold n_{Cells}}' font ",18"
set ylabel '{/:Bold E(t_{End})}' font ",18"
#set yr[0:1.2]
set key bottom left
plot "output/txt/conv/toro3error.txt" using 2:3 title "no AMR (3)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#FF0000" lw 2 lt 2,\
"output/txt/conv/toro3error.txt" using 2:8 title "AMR (3)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#000000" lw 2 lt 2,\
"output/txt/conv/toro3error.txt" using 2:6 title "{/Symbol a}=1" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#0D996A"  lw 2 lt 2,\
"output/txt/conv/toro3error.txt" using 2:7 title "{/Symbol a}=2" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#3355FF"  lw 2 lt 2

set xlabel '{/:Bold n_{Cells}}' font ",18"
set ylabel '{/:Bold E(t_{End})}' font ",18"
#set yr[0:1.2]
set key bottom left
plot "output/txt/conv/toro4error.txt" using 2:3 title "no AMR (4)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#FF0000" lw 2 lt 2,\
"output/txt/conv/toro4error.txt" using 2:8 title "AMR (4)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#000000" lw 2 lt 2,\
"output/txt/conv/toro4error.txt" using 2:6 title "{/Symbol a}=1" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#0D996A"  lw 2 lt 2,\
"output/txt/conv/toro4error.txt" using 2:7 title "{/Symbol a}=2" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#3355FF"  lw 2 lt 2

set xlabel '{/:Bold n_{Cells}}' font ",18"
set ylabel '{/:Bold E(t_{End})}' font ",18"
#set yr[0:1.2]
set origin 0.25,0.01
set key bottom left
plot "output/txt/conv/toro5error.txt" using 2:3 title "no AMR (5)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#FF0000" lw 2 lt 2,\
"output/txt/conv/toro5error.txt" using 2:8 title "AMR (5)" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#000000" lw 2 lt 2,\
"output/txt/conv/toro5error.txt" using 2:6 title "{/Symbol a}=1" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#0D996A"  lw 2 lt 2,\
"output/txt/conv/toro5error.txt" using 2:7 title "{/Symbol a}=2" with linespoints pointtype 5 pointsize 0.6 linecolor rgb "#3355FF"  lw 2 lt 2




set xr[*:*]
set yr[*:*]
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid ls 12
set style line 11 lc rgb '#808080' lt 1
set border 3 ls 11
set key box
set key tmargin
set key font ",10
set key samplen 4  spacing 1.25 
#set output "output/txt/convtoro1.eps"
#replot 
#set term x11
 
 
