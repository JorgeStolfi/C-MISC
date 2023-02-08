#! /bin/bash
# Last edited on 2023-02-04 07:13:41 by stolfi

dfile="$1"; shift
name="${dfile%%.*}"

# Raw data plots
raweps="${name}-raw.eps"
gnuplot <<EOF
set term postscript eps color "TimesRoman" 24
set output "${raweps}"
# set xrange [10:]
# set yrange [10:]
set logscale y
set logscale x
plot \
  "${dfile}" using 1:5 title "ct" with linespoints lt 1 lc rgbcolor '#007700' pt 7, \
  "${dfile}" using 1:6 title "sz" with linespoints lt 2 lc rgbcolor '#ff0000' pt 6
EOF

evince "${raweps}" 

# Binned data plots
tmp=/tmp/$$
bfile="${tmp}-bin.act"
bineps="${name}-bin.eps"

cat ${dfile} \
  | condense-comp-sizes.gawk \
  > ${bfile}

gnuplot <<EOF
set term postscript eps color "TimesRoman" 24
set output "${bineps}"
# set xrange [10:]
# set yrange [10:]
set logscale y
set logscale x
plot \
  "${bfile}" using 2:3 title "ct" with linespoints lt 1 lc rgbcolor '#007700' pt 7, \
  "${bfile}" using 2:4 title "sz" with linespoints lt 2 lc rgbcolor '#ff0000' pt 6
EOF

evince "${bineps}" 

rm -f ${bfile}
