# for reading .csv file

set datafile separator ","

# remove legend
unset key
# set key autotitle columnhead
# plot "< awk '(NR>2){print;}' datafile"

## output size and type
# set terminal epslatex color 10 size 5.2in,2.4in standalone
set terminal pdfcairo color enhanced size 5.2in,2.4in font ", 12pt"
set size ratio 1
set xlabel "p_{11}"
set ylabel "p_{22}"
# set xlabel "$p_{11}$"
# set ylabel "$p_{22}$"

#unset colorbox ## for no colorbox

# # palette
# set palette viridis
set palette rgb 33,13,10

# Plot to 1
set xrange [0.5:1]
show xrange
set yrange [0.5:1]
show yrange
set xtics autofreq 0.5, 0.1, 1 out nomirror font ", 12pt"
set ytics autofreq 0.5, 0.1, 1 out nomirror font ", 12pt"

set output strftime('equilibrium-%H-%M-%S.pdf', time(0))
set multiplot layout 1,2
set title "X_1^*"
plot 'eq-03-19-16-58-19.csv' skip 1 using 1:2:3 w image
set title "X_2^*"
plot 'eq-03-19-16-58-19.csv' skip 1 using 1:2:4 w image
unset multiplot
