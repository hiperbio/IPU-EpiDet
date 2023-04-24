# script adapted from: https://www.xmodulo.com/plot-bar-graph-gnuplot.html

set terminal png size 800,500 enhanced font "Helvetica,20"
set output 'output_chart.png'

red = "#FF0000"; green = "#00FF00"; blue = "#0000FF"; skyblue = "#87CEEB";
set yrange [0:70]
set style data histogram
set style histogram cluster gap 1
set style fill solid
set boxwidth 0.9
set xtics format ""
set grid ytics

set title "M = 6048"
plot "bar.dat" using 2:xtic(1) title "N=4000" linecolor rgb red,	\
     "bar.dat" using 3:xtic(1) title "N=8000" linecolor rgb blue	\

