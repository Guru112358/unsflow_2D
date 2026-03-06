set datafile separator ","
set grid
set xlabel "X"
set ylabel "Cp"
set yrange reverse
set size ratio 1.0
set key outside

set title "RAE 2822, M = 0.725 , AoA= 2.31 degrees"
plot "cp.csv" using 1:2 title "my code" with lines lt rgb "blue" ,"pres.csv" using 1:2 with lines  lt rgb 'red' title "SU2 " 

