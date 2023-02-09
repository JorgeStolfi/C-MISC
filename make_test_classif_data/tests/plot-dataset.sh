#! /bin/bash
# Last edited on 2010-05-24 05:06:59 by stolfi

name="$1"; shift;
nclasses="$1"; shift

cmds="/tmp/$$.gpl"

cat > ${cmds} <<EOF 
  set term postscript eps color linewidth 2 "TimesRoman" 32
  set output "out/${name}-plot.eps"
  set size 3,3
  set size ratio -1
  set title "${name}"
  set xrange [-1.2:+1.2]
  set yrange [-1.2:+1.2]
EOF

printf "  plot" >> ${cmds}
sep=""

ltypes=( 0 1 3 4 5 6 )

class=1
while [[ ${class} -le ${nclasses} ]]; do
  smps=`printf 'out/%s-c%03d.dat' "${name}" ${class}`
  printf '%s \\'"\n" "${sep}" >> ${cmds}
  printf '    "< egrep -e '"'"'^[ ]*[0-9]'"'"' %s" using 2:3 notitle with points lt %d pt 7 pointsize 1.5' ${smps} ${ltypes[${class}]} >> ${cmds}
  sep=","
  class=$(( ${class} + 1 ))
done

printf "\n" >> ${cmds}

cat ${cmds} 2>&1

gnuplot < ${cmds}

gv out/${name}-plot.eps
rm -f ${cmds}
