set terminal pngcairo enhanced

set output "graph1.png"

set log x
set log y

set xlabel "{/Symbol D}t/{/Symbol t}"
set ylabel "L_2-norm"

set title "Temporal error for Taylor-Green vortex on a 64×64 grid"

set key bottom right

plot "fvm_bardow_64.txt" u 1:2 w lp lw 2 t "Bardow", \
    "fvm_dugks_64.txt" u 1:2 w lp lw 2 t "DUGKS", \
     0.001*x w l dt 2 lc rgb "black" t "\\~ {/Symbol D}t"