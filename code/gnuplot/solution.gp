# for reading .csv file

set datafile separator ", "

# remove legend
unset key
# set key autotitle columnhead
# plot "< awk '(NR>2){print;}' datafile"

## output size and type
# set terminal epslatex color 10 size 5.2in,2.4in standalone
set terminal pdfcairo color enhanced size 5.2in,2.4in font "Times, 12pt"
#set size ratio 1
set xlabel "time, [days]"

#unset colorbox ## for no colorbox

# # palette
# set palette viridis
# set palette rgb 33,13,10

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

set xrange [0:730]
#set yrange #[] writeback
#show yrange
set y2tics
#set y2range #[] writeback

set xtics autofreq out nomirror font ", 12pt"
set ytics autofreq out nomirror font ", 12pt"
#et y2tics

set encoding utf8
set output "solution.pdf"

# stats 'patch-dynamics-disjoint-1-1.csv' u "x1" nooutput
# ymax1 = STATS_max
# ymin1 = STATS_min
# ymin1=0.95*ymin1
# ymax1=ymax1*1.05
ymin1=0
ymax1=0.15

set multiplot layout 1,2

set ylabel "X_1"
I1bar = 0.1
# set arrow 1 from 0,I1bar to GPVAL_DATA_X_MAX, I1bar nohead
set arrow 1 from 0,I1bar to 730, I1bar nohead
#set label 1 at 20, I1bar+0.01 '~I_1'
set yrange [ymin1:ymax1]
set y2range [ymin1:ymax1]
set y2tics add ("~I‾_{ 1}" I1bar)
plot 'patch-dynamics-disjoint-1-1.csv' skip 1 using 1:2 w l dt 1 lw 2 lc rgb "red", 'patch-dynamics-irreducible-1.csv' skip 1 using 1:2 w l dt 1 lw 2 lc rgb "blue", 'patch-dynamics-disjoint-1-2.csv' skip 1 using 1:2 w l dt 2 lw 2 lc rgb "red", 'patch-dynamics-irreducible-2.csv' skip 1 using 1:2 w l dt 2 lw 2 lc rgb "blue"
#unset label 1
unset arrow 1

unset yrange
unset y2tics

# stats 'patch-dynamics-disjoint-2-1.csv' u "x2" nooutput
# ymax2 = STATS_max
# ymin2 = STATS_min
# ymin2=0.95*ymin2
# ymax2=ymax2*1.05
ymin2=0
ymax2=0.25

set ylabel "X_2"
I2bar = 0.14
#set arrow 2 from 0,I2bar to GPVAL_DATA_X_MAX, I2bar nohead
set arrow 2 from 0,I2bar to 730, I2bar nohead
#set label 2 at (0, I2bar I2^*)
set yrange [ymin2:ymax2]
set y2range [ymin2:ymax2]
set y2tics add ("~I‾_{ 2}" I2bar)
plot 'patch-dynamics-disjoint-2-1.csv' skip 1 using 1:2 w l dt 1 lw 2 lc rgb "red", 'patch-dynamics-irreducible-1.csv' skip 1 using 1:3 w l dt 1 lw 2 lc rgb "blue", 'patch-dynamics-disjoint-2-2.csv' skip 1 using 1:2 w l dt 2 lw 2 lc rgb "red", 'patch-dynamics-irreducible-2.csv' skip 1 using 1:3 w l dt 2 lw 2 lc rgb "blue"
#unset label 2
unset arrow 2
unset multiplot
