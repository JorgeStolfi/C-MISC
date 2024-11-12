#! /bin/bash
# Last edited on 2023-06-01 10:31:01 by stolfi

# Plots statistics of population size as function of years.

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
yearRange=( `./find_column_range.gawk ${mfile} ${ffile}` )
yStart=${yearRange[0]}
yStop=${yearRange[1]}

# Determine tic spacing for Y axis:
popRange=( `./find_column_range.gawk -v col=4 ${mfile} ${ffile}` )
# echo "popRange = ${popRange[*]}" 1>&2

popMax=${popRange[1]}

echo "popMax = ${popMax}" 1>&2

if [[ ${popMax} -lt 10 ]]; then
  popTics=1
elif [[ ${popMax} -lt 100 ]]; then
  popTics=10
elif [[ ${popMax} -lt 1000 ]]; then
  popTics=100
elif [[ ${popMax} -lt 10000 ]]; then
  popTics=1000
elif [[ ${popMax} -lt 100000 ]]; then
  popTics=10000
else
  echo "** popMax = ${popMax}" 1>&2; exit 1
fi

echo "popTics = ${popTics}" 1>&2

if [[ ${popMax} -lt $(( 2 * ${popTics} )) ]]; then
  popTics=$(( ${popTics} / 10 )); popMTics=2
elif [[ ${popMax} -lt $(( 5 * ${popTics} )) ]]; then
  popTics=$(( ${popTics} / 5 )); popMTics=2
else
  popTics=$(( ${popTics} / 2 )); popMTics=1
fi

echo "popTics = ${popTics} popMTics = ${popMTics}" 1>&2

tmpPlotFile="${tmp}-t.png"
  
export GDFONTPATH="."

gnuplot <<EOF
  set term png truecolor size 1600,1600 font "schlbk,36"
  set output "${tmpPlotFile}"
  set title "max cohort size = ${bpy}"
  
  yStart = ${yStart}
  yStop = ${yStop}
  
  set key right bottom box lc rgb '#ffffff' spacing 1.5 height 0.5
  
  set nologscale y
  set ylabel "population size"
  
  set ytics ${popTics}
  set mytics ${popMTics}
  
  set format y "%.0f"
  set yrange [-0.5:]
  
  load "plot_common.gpl"

  set grid ytics mytics  lt 1 lw 2 lc rgb (grcolor), lt 1 lw 1 lc rgb (grcolor)

  plot \
    "${mfile}" using (s(1)):(v(3)):(v(2)):(v(4)) \
      notitle       with yerrorbars  pt 7 ps 0.0 lw 2 lc rgb '${mcolor}', \
    "${mfile}" using (t(1)):(v(3))\
      title "male"  with linespoints pt 7 ps 1.5 lw 2 lc rgb '${mcolor}', \
    \
    "${ffile}" using (s(1)):(v(3)):(v(2)):(v(4))\
      notitle        with yerrorbars  pt 7 ps 0.0 lw 2 lc rgb '${fcolor}', \
    "${ffile}" using (t(1)):(v(3)) \
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

