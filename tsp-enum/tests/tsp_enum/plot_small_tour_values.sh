#! /bin/bash
# Last edited on 2023-03-31 04:36:07 by stolfi

cmd=${0##*/}
usage="${cmd} -nv NV [-maxsmall MAXSMALL] [-logscale] [-datastyle STYLE] [-show [-timeout NSECS]] FNAME SIDE CVEC..."

# Plots the value-rank plots for some of the cheapest tours, 
# for one or more cost vectors.
# 
# Each {CVEC} should be a distinct single digit, identifying a cost vector.
# The {SIDE} should be "lo" or "hi".
# Creates a single EPS plot file
# showing, for each cost vector {S}, the value-rank plot for the {MAXSMALL}
# cheapest entries in file "{FNAME}-c{S}-asym-{SIDE}.dat", both actual
# and estimated.
#
# The "-logscale" option creates a log-log plot.
# If {MAXSMALL} is not given uses {MAXSMALL=200}
# 

function error()
{ 
  echo "${cmd}: $1" >&2; 
  echo "usage: ${usage}" >&2; 
  exit 1;
}

nv=0
maxsmall=200
logscale=0
datastyle="points";
timeout=0
show=0
while [ $# -gt 0 ]; do
  case "$1" in
    -nv )
      nv=$2; shift; shift ;;
    -maxsmall ) 
      maxsmall=$2; shift; shift ;;
    -logscale ) 
      logscale=1; shift ;;
    -datastyle ) 
      datastyle="$2"; shift; shift ;;
    -show ) 
      show=1; shift ;;
    -timeout ) 
      timeout=$2; shift; shift ;;
    -* )
      error "unrecognized option $1" 
      echo "usage: ${usage}"; exit 1 ;;
    * ) break;;
  esac;
done

if [ ${nv} -le 0 ]; then
  echo 'must specify "-nv"' 
  echo "usage: ${usage}"; exit 1 
fi

if [ $# -lt 3 ]; then
  error "wrong number of parameters"; 
  echo "usage: ${usage}"; exit 1
fi

name="$1"; shift;
side="$1"; shift;
cvecs=( $@ )

tmp="/tmp/$$"

# Create temporary plot files and gather limits:
tpfiles=()

# File with max and min {Z-Zref} values in input files:
zdelfile="${tmp}.zdel"

# File with max and min {T} values in input files:
trangefile="${tmp}.trange"

for cvec in ${cvecs[*]} ; do 
  infile="out/${name}-c${cvec}-asym-${side}.dat"
  tpfile="${tmp}-${cvec}.dat"
  printf "${infile} -> ${tpfile}\n" 1>&2 
  
  # Get the number ${nsmp} of perms in input sample:
  nsmp="`cat ${infile} | egrep -e '^[ ]*[0-9.]' | wc -l | tr -d ' '`"
  printf "input file has ${nsmp} data entries\n" 1>&2 
  
  # Extract the subsample to be plotted:
  cat ${infile} \
    | egrep -e '^[ ]*[0-9.]' \
    | sort -b +2 -3g \
    | head -${maxsmall} \
    > ${tpfile}
  tpfiles=( ${tpfiles[*]} ${tpfile} )

  # Get the number ${nplot} of perms that will be plotted:
  nplot="`cat ${tpfile} | wc -l | tr -d ' '`"
  printf "plotting only ${nplot} cheapest tours\n" 1>&2 
  
  # Find plotting ranges:
  ztmin=(`cat ${tpfile} | head -1`)
  zmin="${ztmin[3]}"
  tmin="${ztmin[2]}"
  ztmax=(`cat ${tpfile} | tail -1`)
  zmax="${ztmax[3]}"
  tmax="${ztmax[2]}"
  printf "observed zmin = %s zmax = %s\n" "${zmin}" "${zmax}" 1>&2 
  printf "expected tmin = %s tmax = %s\n" "${tmin}" "${tmax}" 1>&2 
  printf "%s\n" "${tmin}" >> ${trangefile}
  printf "%s\n" "${tmax}" >> ${trangefile}

  # Extract {Zref} parameter from input file
  zzref=(`cat ${infile} | grep 'Zref' | head -1`)
  zref="${zzref[2]}"
  zdelmin=`gawk -v zr=${zref} -v zm=${zmin} 'BEGIN{printf "%.8f\n", zm - zr;}'`
  printf "%s\n" "${zdelmin}" >> ${zdelfile}
  zdelmax=`gawk -v zr=${zref} -v zm=${zmax} 'BEGIN{printf "%.8f\n", zm - zr;}'`
  printf "%s\n" "${zdelmax}" >> ${zdelfile}
  
  # Remove "e" notation from zref:
  zref=`gawk -v zr=${zref} 'BEGIN{printf "%.8f\n", zr;}'`
  printf "zref = %s zdel = %s _ %s\n" "${zref}" "${zdelmin}" "${zdelmax}" 1>&2 

done

# Get {Z-Zref} range:
zdelmax="`cat ${zdelfile} | sort -b +0 -1gr | head -1`"
zdelmin="`cat ${zdelfile} | sort -b +0 -1g | head -1`"
printf "zdel = %s _ %s\n" "${zdelmin}" "${zdelmax}" 1>&2 

# Get {T} range:
tmax="`cat ${trangefile} | sort -b +0 -1gr | head -1`"
tmin="`cat ${trangefile} | sort -b +0 -1g | head -1`"
printf "t = %s _ %s\n" "${tmin}" "${tmax}" 1>&2 

# Temporary EPS file:
otfile="${tmp}.eps"

# Internal gnuplot script:
gpfile="${tmp}.gnuplot"

# Create gnuplot script:
cat <<EOF > ${gpfile}
set term postscript eps color "TimesRoman" 18
set output
nv=(${nv})
zdelmin=(${zdelmin})
zdelmax=(${zdelmax})
tmin=(${tmin})
tmax=(${tmax})
set title "${fname} ${side}"
set ylabel "T"
set xlabel "Z - Z*"
set key top left
EOF

if [ ${logscale} -ne 0 ]; then
  printf 'set xrange [(zdelmin/1.2):(1.2*zdelmax)]'"\n" >> ${gpfile}
  printf 'set yrange [(tmin/1.5):(tmax*1.5)]'"\n" >> ${gpfile}
  printf 'set logscale'"\n" >> ${gpfile}
else
  printf 'set xrange [(-0.05*zdelmax):(1.15*zdelmax)]'"\n" >> ${gpfile}
  printf 'set yrange [(-0.1*tmax):(1.2*tmax)]'"\n" >> ${gpfile}
  printf 'set yzeroaxis'"\n" >> ${gpfile}
  printf 'set xzeroaxis'"\n" >> ${gpfile}
fi

printf 'plot' >> ${gpfile}
sep=''
for cvec in ${cvecs[*]} ; do 
  infile="${name}-c${cvec}-asym-${side}.dat"
  tpfile="${tmp}-${cvec}.dat"
  printf "issuing plot command for ${tpfile}\n" 1>&2 
  
  # Title for observed data points:
  obstit='title "'"${cvec}"'"'
  
  # Extract {Zref} parameter from input file
  zzref=(`cat ${infile} | grep 'Zref =' | head -1`)
  zref="${zzref[2]}"
  zref=`gawk -v z=${zref} 'BEGIN{printf "%.8f\n",z;}'`
  printf "zref = %s\n" "${zref}" 1>&2 
  
  # Issue command to plot the experimental data points, normalized:
  printf "${sep}"' \\'"\n" >> ${gpfile}
  printf '  "'"${tpfile}"'"' >> ${gpfile} 
  printf ' using (column(4)-('"${zref}"')):(column(3))' >> ${gpfile} 
  printf ' '"${obstit}" >> ${gpfile} 
  printf ' with '"${datastyle}"' linetype '"${cvec}"' pointtype 6' >> ${gpfile}
  
  # Issue command to plot the model-estimated costs and ranks:
  printf   ","' \\'"\n" >> ${gpfile} 
  printf '  "'"${tpfile}"'"' >> ${gpfile} 
  printf ' using (column(6)-('"${zref}"')):(column(5))' >> ${gpfile} 
  printf ' notitle' >> ${gpfile}
  printf ' with lines linetype '"${cvec}" >> ${gpfile}
  sep=','
done

printf "\n" >> ${gpfile}
printf "quit\n" >> ${gpfile}

gnuplot < ${gpfile} \
  | sed -e '/^[%][%]Title[:]/s¦.*¦%%Title: '"${name}"'¦' \
  > ${otfile}

cat ${otfile}

if [ ${show} -ne 0 ]; then
  show-eps -timeout ${timeout} ${otfile}
fi

/bin/rm ${zdelfile} ${trangefile}
/bin/rm ${tpfiles[*]}
/bin/rm ${gpfile}
/bin/rm ${otfile}


