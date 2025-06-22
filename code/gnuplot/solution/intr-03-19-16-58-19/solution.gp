# for reading .csv file

set datafile separator ", "

# remove legend
unset key
# set key autotitle columnhead
# plot "< awk '(NR>2){print;}' datafile"

## output size and type
# set terminal epslatex color 10 size 5.2in,2.4in standalone
#set terminal pdfcairo color enhanced size 5.2in,2.4in font ", 10pt"
set terminal pdfcairo color enhanced size 5.4in,2in font ", 12pt"
set rmargin 8
set size ratio 0.75
set xlabel "t, [дни]"

#unset colorbox ## for no colorbox

# # palette
# set palette viridis
# set palette rgb 33,13,10

## find size of .csv matrix array

# set xtics ("0.5" 0, "0.75" N/2., "1" N-1.)
# set ytics ("0.5" 0, "0.75" M/2., "1" M-1.)

# Places image centered in a box with graduated x and y
# The bos's min and max for x and y are slightly lesser/greater that the true values
# set xrange [] writeback
# set yrange [] writeback

set xrange [0:730]
#set yrange #[] writeback
#show yrange
#set y2range #[] writeback

set xtics autofreq 0, 146, 730 out nomirror font ", 9pt"
set ytics autofreq out nomirror font ", 9pt"
#et y2tics

set encoding utf8
set output strftime('solution-%H-%M-%S.pdf', time(0))

# stats 'patch-dynamics-disjoint-1-1.csv' u "x1" nooutput
# ymax1 = STATS_max
# ymin1 = STATS_min
# ymin1=0.95*ymin1
# ymax1=ymax1*1.05
ymin1=0
ymax1=0.3

set multiplot layout 1,2

set ylabel "X_1(t)"
I1bar = 0.09
# set arrow 1 from 0,I1bar to GPVAL_DATA_X_MAX, I1bar nohead
set arrow 1 from 0,I1bar to 730, I1bar nohead
#set label 1 at 20, I1bar+0.01 '~I_1'
set yrange [ymin1:ymax1]
#set ytics autofreq 0, 146, 730 out nomirror font ", 9pt"
#set ytics auto
#set ytics in nomirror
#set ytics autofreq 0, 146, 730 out nomirror font ", 9pt"
set y2range [ymin1:ymax1]
unset y2tics
set y2tics -1,1.5,2 out nomirror add ("~I‾_{ 1}" I1bar)
#no rep, no mob; no rep, mob; rep, no mob; rep, mob
plot 'patch-dynamics-0.93-0.35-0.85-0.3.csv' skip 1 using 1:2 w l dt 1 lw 1 lc rgb "dark-chartreuse", 'patch-dynamics-0.93-0.35-0.88-0.3.csv' skip 1 using 1:2 w l dt 1 lw 1 lc rgb "midnight-blue", 'patch-dynamics-0.93-0.35-0.92-0.3.csv' skip 1 using 1:2 w l dt 1 lw 1 lc rgb "orchid4", 'patch-dynamics-0.93-0.35-0.97-0.3.csv' skip 1 using 1:2 w l dt 1 lw 1 lc rgb "sienna4"
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
ymax2=0.3

set ylabel "X_2(t)"
I2bar = 0.09
#set arrow 2 from 0,I2bar to GPVAL_DATA_X_MAX, I2bar nohead
set arrow 2 from 0,I2bar to 730, I2bar nohead
#set label 2 at (0, I2bar I2^*)
set yrange [ymin2:ymax2]
set ytics auto
set y2range [ymin2:ymax2]
set y2tics -1,1.5,2 out nomirror add ("~I‾_{ 2}" I2bar)
#no rep, no mob; no rep, mob; rep, no mob; rep, mob
plot 'patch-dynamics-0.93-0.35-0.85-0.3.csv' skip 1 using 1:3 w l dt 1 lw 1 lc rgb "dark-chartreuse", 'patch-dynamics-0.93-0.35-0.88-0.3.csv' skip 1 using 1:3 w l dt 1 lw 1 lc rgb "midnight-blue", 'patch-dynamics-0.93-0.35-0.92-0.3.csv' skip 1 using 1:3 w l dt 1 lw 1 lc rgb "orchid4", 'patch-dynamics-0.93-0.35-0.97-0.3.csv' skip 1 using 1:3 w l dt 1 lw 1 lc rgb "sienna4"
#unset label 2
unset arrow 2
unset multiplot
