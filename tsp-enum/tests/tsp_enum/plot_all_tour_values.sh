#! /bin/bash
# Last edited on 2023-03-31 04:35:55 by stolfi

cmd=${0##*/}
usage="${cmd} -nv NV [-maxplot MAXPLOT] [-datastyle STYLE] [-show [-timeout NSECS]] FNAME CVEC"

# Plots costs for perms on {NV} vertices. 
# The {CVEC} should be a single digit, identifying an arc cost vector. 
# Creates an EPS
# plot file showing the value-rank plot for the entries in file
# "{FNAME}-c{CVEC}.dat".
# 
# If {MAXPLOT} is specified, plots only a subsample of {MAXPLOT} entries,
# including the cheapest and most expensive ones.

function error()
{ 
  echo "${cmd}: $1" >&2; 
  echo "usage: ${usage}" >&2; 
  exit 1;
}

nv=0
maxplot=99999999
datastyle="points";
timeout=0
show=0
while [ $# -gt 0 ]; do
  case "$1" in
    -nv )
      nv=$2; shift; shift ;;
    -maxplot ) 
      maxplot="$2"; shift; shift ;;
    -datastyle ) 
      datastyle="$2"; shift; shift ;;
    -show ) 
      show=1; shift ;;
    -timeout ) 
      timeout="$2"; shift; shift ;;
    -* )
      error "unrecognized option $1"; exit 1 ;;
    * ) break;;
  esac;
done

if [ ${nv} -le 0 ]; then
  echo 'must specify "-nv"' 
  echo "usage: ${usage}"; exit 1 
fi

if [ $# -ne 2 ]; then
  error "wrong number of parameters"; 
fi

name="$1"; shift;
cvec="$1"; shift;

tmp="/tmp/$$"

# Get input data file:
infile="out/${name}-c${cvec}-glob.dat"

# Get the number ${nsmp} of perms in input sample:
nsmp="`cat ${infile} | egrep -e '^[ ]*[0-9.]' | wc -l | tr -d ' '`"
printf "read data of ${nsmp} tours\n" 1>&2 

# Create temporary plot file and gather limits:
tpfile="${tmp}-c${cvec}.tmp"
printf "${infile} -> ${tpfile}\n" 1>&2 

# Decide max number ${npick} of perms to plot:
if [ ${nsmp} -lt ${maxplot} ]; then npick=${nsmp} ; else npick=${maxplot} ; fi

# Extract the subsample to be plotted, and fit approximate formula:
cat ${infile} \
  | egrep -e '^[ ]*[0-9.]' \
  | sort -b +2 -3g \
  | pick-some-records -v nin=${nsmp} -v nout=${npick} \
  > ${tpfile}

# Get the number ${nplot} of perms that will be plotted:
nplot="`cat ${tpfile} | wc -l | tr -d ' '`"
printf "plotting only ${nplot} equally spaced tours\n" 1>&2 

# Find plotting ranges:
ztmin=(`cat ${tpfile} | head -1`)
zmin="${ztmin[3]}"
tmin="${ztmin[4]}"
ztmax=(`cat ${tpfile} | tail -1`)
zmax="${ztmax[3]}"
tmax="${ztmax[4]}"
printf "observed zmin = %s zmax = %s\n" "${zmin}" "${zmax}" 1>&2 
printf "expected tmin = %s tmax = %s\n" "${tmin}" "${tmax}" 1>&2 

# Temporary EPS file:
otfile="${tmp}.eps"

# Internal gnuplot script:
gpfile="${tmp}.gnuplot"

# Create gnuplot script:
cat <<EOF > ${gpfile}
set term postscript eps color "TimesRoman" 18
set output
zmin=(${zmin})
zmax=(${zmax})
zdel=(zmax-zmin)
nplot=(${nplot})
tmin=0.0
tmax=1.0
tdel=(tmax-tmin)
set title "${fname} ${cvec}"
set ylabel "T"
set xlabel "Z"
set key top left
EOF

printf 'set xrange [(zmin-0.05*zdel):(zmax+0.15*zdel)]'"\n" >> ${gpfile}
printf 'set yrange [(tmin-0.1*tdel):(tmax+0.2*tdel)]'"\n" >> ${gpfile}
printf 'set yzeroaxis'"\n" >> ${gpfile}
printf 'set xzeroaxis'"\n" >> ${gpfile}

printf 'plot' >> ${gpfile}
sep=''

printf "issuing plot command for ${tpfile}\n" 1>&2 

# Title for observed data points:
obstit='title "'"c${cvec}"'"'

# Issue command to plot the actual costs and relative ranks:
printf "${sep}"' \\'"\n" >> ${gpfile}
printf '  "'"${tpfile}"'"' >> ${gpfile} 
printf ' using (column(4)):(column(3))' >> ${gpfile} 
printf ' '"${obstit}" >> ${gpfile} 
printf ' with '"${datastyle}"' linetype 1 pointtype 6' >> ${gpfile}

sep=','

# Issue command to plot the estimated costs and relative ranks:
printf "${sep}"' \\'"\n" >> ${gpfile}
printf '  "'"${tpfile}"'"' >> ${gpfile} 
printf ' using (column(6)):(column(5))' >> ${gpfile} 
printf ' notitle' >> ${gpfile}
printf ' with lines linetype 3' >> ${gpfile}

sep=','

# Issue command to plot the upper limit {t(Z)=1}:
printf "${sep}"' \\'"\n" >> ${gpfile}
printf ' (1)' >> ${gpfile} 
printf ' notitle' >> ${gpfile} 
printf ' with lines linetype 0' >> ${gpfile}

sep=','

printf "\n" >> ${gpfile}
printf "quit\n" >> ${gpfile}

gnuplot < ${gpfile} \
  | sed -e '/^[%][%]Title[:]/s¦.*¦%%Title: '"${name}-${net}"'¦' \
  > ${otfile}

cat ${otfile}

if [ ${show} -ne 0 ]; then
  show-eps -timeout ${timeout} ${otfile}
fi

/bin/rm ${tpfile}
/bin/rm ${gpfile}
/bin/rm ${otfile}
