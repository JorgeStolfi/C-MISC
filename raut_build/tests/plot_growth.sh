#! /bin/bash
# Last edited on 2009-10-31 00:22:58 by stolfi

prefix="$1"; shift;
pltfile="${prefix}.plt"
pngfile="${prefix}.png"

rm -f ${pngfile}
gnuplot <<EOF
set term png small size 600,800 \
  xffffff x000000 x404040 \
  xff0000 x008000 x0054ff xa800ff \
  x804000 x006780 x5426ff x930080
set output "${pngfile}"
plot \
  "${pltfile}" using 1:4 title "max" with linespoints, \
  "${pltfile}" using 1:6 title "root" with linespoints, \
  "${pltfile}" using 1:7 title "reach" with linespoints
quit
EOF
if [[ -s ${pngfile} ]]; then
  display ${pngfile}
fi  
