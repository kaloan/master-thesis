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

## find size of .csv matrix array

M = 0   # number of rows, 0 if unknown
N = 0   # number of cols, 0 if unknown

if (!M || !N) {
    stats 'eq2.csv' u 0 nooutput   # get number of rows and cols
    M = STATS_records
    N = STATS_columns
}

#set cbtics autofreq

# Plot to 1
set xrange [0.8:1]
show xrange
set yrange [0.8:1]
show yrange
set xtics autofreq 0.8, 0.05, 1 out nomirror font ", 12pt"
set ytics autofreq 0.8, 0.05, 1 out nomirror font ", 12pt"

set output "equilibrium1.pdf"
set multiplot layout 1,2
set title "X_1^*"
plot 'eq4.csv' skip 1 using 1:2:3 w image
set title "X_2^*"
plot 'eq4.csv' skip 1 using 1:2:4 w image
unset multiplot

# Plot to 0.98
set xrange [0.8:0.98]
show xrange
set yrange [0.8:0.98]
show yrange
set xtics autofreq 0.8, 0.06, 0.98 out nomirror font ", 12pt"
set ytics autofreq 0.8, 0.06, 0.98 out nomirror font ", 12pt"

set output "equilibrium2.pdf"
set multiplot layout 1,2
set cbrange [0:0.0220246810603182]
set title "X_1^*"
plot 'eq4.csv' skip 1 using 1:2:3 w image
set cbrange [0.00274363372857308:0.0697692196950718]
set title "X_2^*"
plot 'eq4.csv' skip 1 using 1:2:4 w image
unset multiplot

# Plot to 0.95
set xrange [0.8:0.95]
show xrange
set yrange [0.8:0.95]
show yrange
set xtics autofreq 0.8, 0.05, 0.95 out nomirror font ", 12pt"
set ytics autofreq 0.8, 0.05, 0.95 out nomirror font ", 12pt"

set output "equilibrium3.pdf"
set multiplot layout 1,2
set cbrange [0:0.0215238894548882]
set title "X_1^*"
plot 'eq4.csv' skip 1 using 1:2:3 w image
set cbrange [0.00274363372857308:0.0566752985695734]
set title "X_2^*"
plot 'eq4.csv' skip 1 using 1:2:4 w image
