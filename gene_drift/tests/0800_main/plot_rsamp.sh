#! /bin/bash
# Last edited on 2023-05-31 21:23:08 by stolfi

# Plots ratios femi/masc of surviving lineages per generation,
# for several runs with one population size

# Arguments are {outFullPrefix} {bpy}
# where {outFullPrefix} if the fileprefix including the TeX tag ("A", "B", ..")
# and {bpy} is the max cohort size.

# Reads the file "{outFullPrefix}-rsamp.txt".
# Assumes that the file is the output of {gdr_main.c}.

outFullPrefix="$1"; shift;  # Output file prefix, with tag, minus "-rsamp.txt".
bpy="$1"; shift    # Population size, for title.           

tmp="/tmp/$$"

rfile="${outFullPrefix}-rsamp.txt"
rcolor="#338833"

# Determine first and last years:
yRange=( `./find_column_range.gawk ${rfile}` )
echo "yRange = ${yRange[*]}" 1>&2
yStart=${yRange[0]}
yStop=${yRange[1]}

tmpPlotFile="${tmp}-t.png"
 
export GDFONTPATH="."

gnuplot <<EOF
  set term png truecolor size 2400,1600 font "schlbk,36"
  set output "${tmpPlotFile}"
  set title "female/male ratios - max cohort size = ${bpy}"
  
  yStart = ${yStart}
  yStop = ${yStop}

  set nokey
  
  set nologscale y
  set yrange [-0.01:30.01]
  set ylabel "female/male surviving lineage ratio rsl(g)"
  
  set ytics 2
  set mytics 2
  
  load "plot_common.gpl"
  
  set grid ytics mytics  lt 1 lw 2 lc rgb (grcolor), lt 1 lw 1 lc rgb (grcolor)

  ncolors = 10
  rgb(R,G,B) = 65536 * int(R*255) + 256 * int(G*255) + int(B*255)
  hue(z) = rgb(0.60 + 0.30*z, 0.50 - 0.20*z, 0.20)
  clr(kcol) = hue((int(column(kcol)) % ncolors + 0.0)/(ncolors - 1.0))
  
  plot \
    "${rfile}" using (t(1)):3:(clr(2)) notitle with linespoints pt 7 ps 0.7 lw 2.0 lc variable
  quit
EOF

plotFile="${outFullPrefix}-rsamp.png"

if [[ -s ${tmpPlotFile} ]]; then
  convert ${tmpPlotFile} -resize '50%' ${plotFile}
  display ${plotFile}
else
  echo "** plot file not generated" 1>&2 ; exit 1
fi

rm -f ${tmp}*.png

