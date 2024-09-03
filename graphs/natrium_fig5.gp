reset
set terminal pngcairo enhanced font ",12"
set output "fig5_logx.png"


set xtics format "%3.0e"
set xlabel "Time Step dt"
set ylabel "Error <|u - uref|_2>/U

set xrange [1.0e-8:1.0e-6]
set yrange [0.00007:0.00021]

set log x

z1(x) = x > 5.0e-7 ? 140.0*x : 1/0
z2(x) = x > 5.0e-7 ? 2.0e8*x**2 : 1/0

file = "natrium_fig5.txt"

plot file u 1:2 w lp pt 6 t "SLLBM, p=4",\
z1(x) t "O(dt)", \
z2(x) t "O(dt^2)"

# ----

set output "fig5_loglog.png"
set log y

set ytics add ("0.00020" 0.0002)

replot