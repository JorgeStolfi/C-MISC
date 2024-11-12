#! /bin/bash 
# Last edited on 2011-08-03 23:57:36 by stolfilocal

pop="$1"; shift
gen="$1"; shift
dat_file="$1"; shift

tmp="/tmp/$$"

tmp_plot="${tmp}-tmp.png"
red_plot="${tmp}-red.png"

echo "plotting..." 1>&2

set GDFONTPATH="."

gnuplot << EOF
  set term png size 1600,1200 font "arial" 18
  set output "${tmp_plot}" 
  set lmargin 8; set rmargin 4; set tmargin 0.5; set bmargin 4.0
  set xlabel "Generation"
  set ylabel "Number of RED individuals"
  unset title
  # set grid ytics nomytics
  set bars 0.0
  pop = ${pop}
  gen = ${gen}
  set yrange [-1:(pop+1)]; 
  set xrange [-1:(gen+1)]; 
  # col(k) = (column(0)==0 ? column(k) : 0/0)
  col(k) = column(k)
  plot "${dat_file}" using 2:(col(3)) notitle with lines lc rgb '#00aa00'
EOF

echo "reducing..." 1>&2

convert ${tmp_plot} -resize '50%' ${red_plot}

display ${red_plot}

rm -f ${tmp}-*

echo "done." 1>&2
