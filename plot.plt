set datafile separator ","
set title "Convergence history"
set grid
set mouse
set xlabel "Iterations"
set ylabel "Residual"
set logscale xy
set autoscale xy
plot "convergence_history.csv" using 1:2 with lines  title "rho",\
       "convergence_history.csv" using 1:3 with lines title "rhoU",\
       "convergence_history.csv" using 1:4 with lines title  "rhoV",\
       "convergence_history.csv" using 1:5 with lines title  "rhoE"
       

pause 20
reread
