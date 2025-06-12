# for reading .csv file

set datafile separator ","

# remove legend
unset key
# set key autotitle columnhead
# plot "< awk '(NR>2){print;}' datafile"

## output size and type
set terminal epslatex color 10 size 5.2in,2.4in standalone
# set terminal pdfcairo color enhanced size 5.2in,2.4in
set size ratio 1
# set xlabel "p_{11}"
# set ylabel "p_{22}"
set xlabel "$p_{11}$"
set ylabel "$p_{22}$"

#unset colorbox ## for no colorbox

# # palette
set palette viridis

## find size of .csv matrix array

M = 0   # number of rows, 0 if unknown
N = 0   # number of cols, 0 if unknown

if (!M || !N) {
    stats 'eq2.csv' u 0 nooutput   # get number of rows and cols
    M = STATS_records
    N = STATS_columns
}


# set xtics ("0.5" 0, "0.75" N/2., "1" N-1.)
# set ytics ("0.5" 0, "0.75" M/2., "1" M-1.)

# Places image centered in a box with graduated x and y
# The bos's min and max for x and y are slightly lesser/greater that the true values
# set xrange [] writeback
# set yrange [] writeback

set xrange [0.8:1]
show xrange
set yrange [0.8:1]
show yrange
set tics out nomirror font ", 12"

set output "equilibrium.tex"
# set output "equilibrium.pdf"
set multiplot layout 1,2
set title "$X_1^*$"
plot 'eq4.csv' skip 1 using 1:2:3 w image
set title "$X_2^*$"
plot 'eq4.csv' skip 1 using 1:2:4 w image
# plot 'eq.csv' skip 1 matrix u 1:2:3 w image
# plot 'eq.csv' skip 1 matrix u 1:2:4 w image

# unset output

# # for parameters for I1max, I2max
# I1max = 25
# I2max = 45

# # # 2 plots, I1* < I1max, I2* < I2max:
# # multiplot
# #set output "ex2.tex"
# set output "equilibrium2.pdf"
# set multiplot layout 1,2

# set title "$I_1^{\\ast} \\le \\hat I_1$"
# set title "I_1^* <I_1"
# plot 'exampleplot.csv' matrix u 1:2:($3<I1max) w image palette

# set title "$I_2^{\\ast} \\le \\hat I_2$"
# plot 'exampleplot.csv' matrix u 1:2:($3<I2max) w image palette

# unset multiplot

# unset output

# # # 3 plots, I1* < I1max, I2* < I2max, I1* < I1max & I2* < I2max

# array mtrx[M*N]
# stats 'exampleplot.csv' matrix u (mtrx[int($2*N+$1+1)] = ($3<I1max) ) nooutput

# set output "ex3.tex"#
# set output "ex3.pdf"
# set multiplot layout 1,3

# set title "$I_1^{\\ast} \\le \\hat I_1$"
# plot 'exampleplot.csv' matrix u 1:2:($3<I1max) w image

# set title "$I_2^{\\ast} \\le \\hat I_2$"
# plot 'exampleplot.csv' matrix u 1:2:($3<I2max) w image

# set title "$I_1^{\\ast} \\le \\hat I_1\\wedge I_2^{\\ast} \\le \\hat I_2$"
# plot 'exampleplot.csv' matrix u 1:2:(($3<I2max)&&mtrx[int($2*N+$1+1)]) w image

# unset multiplot

# unset output
