#! /bin/bash
# Last edited on 2012-12-20 22:23:28 by stolfilocal

datafile="$1"; shift;

name="${datafile%%.*}"
echo "${name}" 1>&2

gnuplot <<EOF
set term x11
set yrange[1.0e-10:]
set title "${name##*/}"
set logscale y
plot \
  "${datafile}" using 2:6 title "fmin" with linespoints, \
  "" using 2:7 title "fmax" with linespoints, \
  "" using 2:8 title "rdif" with linespoints
pause 120
EOF
