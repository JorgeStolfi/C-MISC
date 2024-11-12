#! /bin/bash
# Last edited on 2023-05-31 21:22:51 by stolfi

# Plots ratio femi/masc of surviving lineages per generation,
# for two values of the max cohort size {bpy}.

# Arguments are {outPrefix} {bpy[1]}.{tag[1]} {bpy[2].tag[2]} ...
# where {bpy[r]} is the max cohort size of program run {r},
# and {tag[r]}is its TeX tag.

# Reads the files "{outPrefix}-{texTag}-ratio.txt".
# Assumes that the file is the output of {gdr_main.c}.

outPrefix="$1"; shift;  # Output file prefix, minus "-{f,m}-nlins.txt".

bpys=()
tags=()
while [[ $# -ge 1 ]]; do
  pt="$1"; shift
  bpys+=( "${pt/.*/}" )
  tags+=( "${pt/*./}" )
done

if [[ ${#bpys[@]} -ne 2 ]]; then
  echo "** can handle only 2 plots for now" 1>&2; exit 1
fi

tmp="/tmp/$$"

rfilea="${outPrefix}-${tags[0]}-ratio.txt"
rcolora="#008800"

rfileb="${outPrefix}-${tags[1]}-ratio.txt"
rcolorb="#8800ff"

# Determine first and last years:
yRange=( `./find_column_range.gawk ${rfilea} ${rfileb}` )
yStart=${yRange[0]}
yStop=${yRange[1]}

tmpPlotFile="${tmp}-t.png"
 
export GDFONTPATH="."

gnuplot <<EOF
  set term png truecolor size 2400,1600 font "schlbk,36"
  set output "${tmpPlotFile}"
  set title "female/male ratios"
  
  yStart = ${yStart}
  yStop = ${yStop}

  set key left top reverse Left box lc rgb '#ffffff' spacing 1.5 height 0.5
  
  set nologscale y
  set yrange [-0.01:30.01]
  set ylabel "median and range of female/male lineage ratio rsl(g)"
  
  set ytics 2
  set mytics 2
  
  load "plot_common.gpl"
  
  set grid ytics mytics  lt 1 lw 2 lc rgb (grcolor), lt 1 lw 1 lc rgb (grcolor)

  plot \
    "${rfilea}" using (s(1)):3:2:4  notitle                  with yerrorbars  pt 7 ps 0.0 lw 2.5 lc rgb '${rcolora}', \
    "${rfilea}" using (t(1)):3      title "bpy = ${bpys[0]}" with linespoints pt 7 ps 1.5 lw 2.0 lc rgb '${rcolora}', \
    \
    "${rfileb}" using (s(1)):3:2:4  notitle                  with yerrorbars  pt 7 ps 0.0 lw 5.0 lc rgb '${rcolorb}', \
    "${rfileb}" using (t(1)):3      title "bpy = ${bpys[1]}" with linespoints pt 7 ps 1.5 lw 2.0 lc rgb '${rcolorb}'
  quit
EOF

plotFile="${outPrefix}-ratios.png"

if [[ -s ${tmpPlotFile} ]]; then
  convert ${tmpPlotFile} -resize '50%' ${plotFile}
  display ${plotFile}
else
  echo "** plot file not generated" 1>&2 ; exit 1
fi

rm -f ${tmp}*.png

