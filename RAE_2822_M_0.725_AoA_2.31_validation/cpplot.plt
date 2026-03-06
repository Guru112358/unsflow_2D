set datafile separator ","
set grid
set mouse
set xlabel "X"
set ylabel "Cp"
set yrange reverse
set size ratio 1
set key outside
set autoscale xy
set title "RAE 2822, M = 0.85 , AoA= 2.31 degrees"
plot "cp.csv" using 1:2 title "my code" with lines lt rgb "blue" ,"pres.csv" using 1:2 with lines  lt rgb 'red' title "SU2 benchmark" 

