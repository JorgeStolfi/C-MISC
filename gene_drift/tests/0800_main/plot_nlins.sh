#! /bin/bash
# Last edited on 2023-06-01 10:31:11 by stolfi

# Plots number of surviving lineages per generation, "m" and "f",
# as percentages of the correponding population.

# Reads the files "{fullPrefix}-{sex}-nlins.txt" where {sex} is "0" or "1".
# Assumes that each file is the output of {gdr_main.c}.

fullPrefix="$1"; shift; # Output file prefix, with the TeX tag, minus "-{0,1}-nlins.txt".
bpy="$1"; shift;        # Max births per year (for title).

tmp="/tmp/$$"

mfile="${fullPrefix}-0-nlins.txt"
mcolor="#0022ff"

ffile="${fullPrefix}-1-nlins.txt"
fcolor="#ff2200"

# Determine first and last years:
yRange=( `./find_column_range.gawk ${mfile} ${ffile}` )
yStart=${yRange[0]}
yStop=${yRange[1]}

tmpPlotFile="${tmp}-t.png"
  
export GDFONTPATH="."

gnuplot <<EOF
  set term png truecolor size 1600,1600 font "schlbk,36"
  set output "${tmpPlotFile}"
  set title "max cohort size = ${bpy}"
  
  yStart = ${yStart}
  yStop = ${yStop}
  
  set key left bottom reverse Left box lc rgb '#ffffff' spacing 1.5 height 0.5
  
  set logscale y
  set ylabel "surviving lineages (% of population)"
  
  set ytics ( \
    0.001 0, 0.002 0, "" 0.003 0, "" 0.004 0, 0.005 0, "" 0.006 0, "" 0.007 0, "" 0.008 0, "" 0.009 0, \
    0.01 0, 0.02 0, "" 0.03 0, "" 0.04 0, 0.05 0, "" 0.06 0, "" 0.07 0, "" 0.08 0, "" 0.09 0, \
    0.1 0, 0.2 0, "" 0.3 0, "" 0.4 0, 0.5 0, "" 0.6 0, "" 0.7 0, "" 0.8 0, "" 0.9 0, \
    1 0, 2 0, "" 3 0, "" 4 0, 5 0, "" 6 0, "" 7 0, "" 8 0, "" 9 0, \
    10 0, 20 0, "" 30 0, "" 40 0, 50 0, "" 60 0, "" 70 0, "" 80 0, "" 90 0, \
    100 0, 200 0, "" 300 0, "" 400 0, 500 0, "" 600 0, "" 700 0, "" 800 0, "" 900 0 \
  )
  set format y "%.2f"
  set yrange [0.009:101.0]
  
  load "plot_common.gpl"

  set grid ytics mytics  lt 1 lw 2 lc rgb (grcolor), lt 1 lw 1 lc rgb (grcolor)

  plot \
    "${mfile}" using (s(1)):(v(6)):(v(5)):(v(7)) \
      notitle       with yerrorbars  pt 7 ps 0.0 lw 2 lc rgb '${mcolor}', \
    "${mfile}" using (t(1)):(v(6))\
      title "male"  with linespoints pt 7 ps 1.5 lw 2 lc rgb '${mcolor}', \
    \
    "${ffile}" using (s(1)):(v(6)):(v(5)):(v(7))\
      notitle        with yerrorbars  pt 7 ps 0.0 lw 2 lc rgb '${fcolor}', \
    "${ffile}" using (t(1)):(v(6)) \
      title "female" with linespoints pt 7 ps 1.5 lw 2 lc rgb '${fcolor}'
  quit
EOF

plotFile="${fullPrefix}-nlins.png"

if [[ -s ${tmpPlotFile} ]]; then
  convert ${tmpPlotFile} -resize '50%' ${plotFile}
  display ${plotFile}
else
  echo "** plot file not generated" 1>&2 ; exit 1
fi

rm -f ${tmp}*.png

