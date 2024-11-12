#! /bin/bash
# Last edited on 2020-01-08 13:00:03 by jstolfi

TAG="$1"; shift   # Tag for output files.
N="$1"; shift     # Number of terms on each side of center one.

logy=0
show=YES

export GDFONTPATH=../ttf

temp="/tmp/$$"
dfile="out/funcs_${TAG}.txt"
ffile="${temp}-full.png"
gfile="${temp}.gpl"
pfile="out/funcs_${TAG}.png"

if [[ ! ( -s ${dfile} ) ]]; then
  echo "** ${dfile} missing or empty" 1>&2; exit 1
fi

# Build the plot command:
printf "plot \\\\\n" > ${gfile}
printf "'${dfile}' using 2:3        title 'wt'    with lines lt 1 lc rgb '#666666', \\\\\n" >> ${gfile}
printf "'${dfile}' using 2:4        title 'H'     with lines lt 1 lc rgb '#ff0088', \\\\\n" >> ${gfile}
nxcol=5 # Next data column.
kk=$(( 0 - ${N} )) # Term index, {-N..+N}.
while [[ ${kk} -le ${N} ]] ; do 
  if [[  kk -lt 0 ]]; then k=$(( -${kk} )) ; else k=${kk} ; fi # Coeff index,{0..N}.
  printf "'${dfile}' using 2:${nxcol} title 'C${k}*G[${kk}]' with lines lt 1 lc rgb '#008800', \\\\\n" >> ${gfile}
  nxcol=$(( ${nxcol} + 1 ))
  kk=$(( ${kk} + 1 ))
done
fcol=${nxcol} # Data column of {F} values.
ecol=$(( ${fcol} + 1 )) # Data column of {F-H} values.
printf "'${dfile}' using 2:${fcol}        title 'F'     with lines lt 1 lc rgb '#0055ff', \\\\\n" >> ${gfile}
printf "'${dfile}' using 2:(err(${ecol})) title 'F-H'   with lines lt 1 lc rgb '#ff7700'\n"  >> ${gfile}


gnuplot <<EOF
set term png size 2400,1200 font "courbd,24"
set output "${ffile}"

set xtics mirror 1.0
set mxtics 10

if (${logy} > 0) {
  set logscale y
  set yrange [-0.9e-12: 1.1]
  err(k) = abs(column(k))
} 
else {
  set yrange [-0.61 : 0.61 ]
  set ytics mirror 0.1
  set mytics 10
  set yzeroaxis lt 1 lw 3 lc rgb '#ffddaa'
  err(k) = column(k)
}
set grid xtics  lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 2 lc rgb '#ffddaa'
set grid ytics  lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 2 lc rgb '#ffddaa'

set xzeroaxis lt 1 lw 3 lc rgb '#ffddaa'

load "${gfile}"
EOF

if [[ -s ${ffile} ]]; then 
  convert ${ffile} -resize '50%' ${pfile}
  if [[ "/${show}" == "/YES" ]]; then display ${pfile}; fi
fi

rm -fv ${temp}{-*,}.*
