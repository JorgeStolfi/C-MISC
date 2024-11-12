#! /bin/bash
# Last edited on 2015-04-21 22:14:59 by stolfilocal

# Plots the prices of one or more exchanges.
# 
# Arguments: {TITLE} {TTICDAYS} {TMTICS} {DMIN} {DMAX} {VMIN} {VMAX} {LOGSCALE} {SHOW} \
#    {FNAME[0]}   {FIX[0]}   {LABEL[0]}   {RATE[0]}   {COLOR[0]} \
#    {FNAME[1]}   {FIX[1]}   {LABEL[1]}   {RATE[1]}   {COLOR[1]} \
#    ... \
#    {FNAME[N-1]} {FIX[N-1]} {LABEL[N-1]} {RATE[N-1]} {COLOR[N-1]} 
# 
# Where 
#   
#   {TITLE} is the plot's title;
#   {TTICDAYS} is the major tic spacing on the X-axis, in days;
#   {TMTICS} is the number of sub-intervals per major time tic interval, defined by minor tics and grid; 
#   {DMIN,DMAX} is the range of dates to plot (format "%Y-%m-%d");
#   {VMIN},{VMAX} is the range of prices for the vertical scale of the plot (or both 0 to leave it free); 
#   {LOGSCALE} is 0 for linear scale, 1 for logscale;
#   {SHOW} is "YES" to display the plot.
# 
#   {FNAME[i]} is a data file with summary trade data (open, high, low etc.); 
#   {FIX[i]} is the index of the field in {FNAME[i]} to be plotted (from 1); 
#   {LABEL[i]} is the label to show on the plot's key; 
#   {RATE[i]} is a factor to be divided into the price as read from the file.
#   {COLOR[i]} is the RGB color to be used for this file (6 hex digits).
#   
# Each input file "{FNAME[i]}" must contain one line per time interval.
# The first two fields must be "{DATE} {TIME}". 
#
# If there is only one file, the key is suppressed, so be sure to pick the right title.

title="$1"; shift  # Plot title.
tticdays="$1"; shift  # Number of days per tic interval (may be fractional).
tmtics="$1"; shift    # Number of minor tic intervals per major tic interval.
dmin="$1"; shift      # Min date to plot.
dmax="$1"; shift      # Max date to plot.
vmin="$1"; shift      # Min price for plot range.
vmax="$1"; shift      # Max price for plot range.
logScale="$1"; shift  # Linear scale (0) or log scale (1).
show="$1"; shift      # "YES" to display plot.
fname=()
fix=()
label=()
rate=()
color=()
i=0;
while [[ $# -gt 1 ]]; do
  fname[$i]="$1"; shift # File name.
  fix[$i]="$1"; shift   # Field index.
  label[$i]="$1"; shift # Plot label.
  rate[$i]="$1"; shift  # Currency conversion factor.
  color[$i]="$1"; shift # Color for plotting.
  i=$(( $i + 1 ))
done
nfiles=$i

tmp="/tmp/$$"

gplFile="${tmp}.gpl"

fullplotFile="${tmp}-full.png"
plotFile="${tmp}.png"

# Create the plot command
printf "plot" > ${gplFile}
sep=""

echo "plotting files:" 1>&2
i=0;
while [[ $i -lt $nfiles ]]; do
  fnamei="${fname[$i]}"
  fixi=${fix[$i]}
  ratei=${rate[$i]}
  colori="${color[$i]}"
  
  if [[ "/${ratei}" =~ /[1]([.][0]*)? ]]; then 
    titlei="${label[$i]}"
  else
    titlei="${label[$i]}/${ratei}"
  fi
  echo "  column ${fixi} of ${fnamei}" 1>&2
  printf "%s "'\\'"\n" "${sep}" >> ${gplFile}
  printf "  \"${fnamei}\" using (tim(1)):(pmd(${fixi},${ratei})) title \"${titlei}\"" >> ${gplFile}
  printf " with linespoints pt 7 ps 0.75 lt 1 lw 1.5 lc rgb '#${colori}'" >> ${gplFile}
  sep=","
  i=$(( $i + 1 ))
done
printf "\n" >> ${gplFile}

# cat ${gplFile} 1>&2 

export GDFONTPATH=.:..

gnuplot <<EOF
set term png size 3000,1500 font "courbd,24"
set output "${fullplotFile}"
logScale=${logScale}
nfiles=${nfiles}

tticdays = ${tticdays}       # Number of days per time tic interval
tticsecs = 60*60*24*tticdays # Number of seconds per time tic interval
tmtics =  ${tmtics}          # Number of minor time tic intervals per major time tic interval

set xdata time
xmin=strptime("%Y-%m-%d","${dmin}")
xmax=strptime("%Y-%m-%d","${dmax}")
set xrange [(xmin-24*3600):(xmax+24*3600)]

unset xtics
set xtics mirror tticsecs
set format x "%Y-%m-%d"
set mxtics tmtics

set title "${title}"

# Input time:
set timefmt "%Y-%m-%d %H:%M:%S"
tim(kcol) = (timecolumn(kcol))

# Price field, adjusted by given rate:
pmd(kpmd,rate) = (column(kpmd) == 0 ? 0/0 : column(kpmd)/rate)

vmin = ${vmin}
vmax = ${vmax}

if (vmin < vmax) {
  set yrange [vmin:vmax]
}
  
if (logScale > 0) {
  set logscale y

  if (vmax/vmin > 100) {
    set ytics ( \
      0.010,0.020,0.030,0.050,0.070, \
      0.10, 0.20, 0.30, 0.50, 0.70,  \
      1,    2,    3,    5,    7,     \
      10,   20,   30,   50,   70,    \
      100,  200,  300,  500,  700,   \
      1000, 2000, 3000, 5000, 7000   \
    )
  } else {
  if (vmax/vmin > 10) {
    set ytics ( \
      0.010,0.012,0.015,0.020,0.025,0.030,0.040,0.050,0.060,0.070,0.080,0.090, \
      0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90,  \
      1,    1.2,  1.5,  2,    2.5,  3,    4,    5,    6,    7,    8,    9,     \
      10,   12,   15,   20,   25,   30,   40,   50,   60,   70,   80,   90,    \
      100,  120,  150,  200,  250,  300,  400,  500,  600,  700,  800,  900,   \
      1000, 1200, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 7000, 8000, 9000   \
    )
  } else {
    set ytics ( \
      0.010,0.011,0.012,0.015,0.020,0.025,0.030,0.035,0.040,0.045,0.050,0.055,0.060,0.065,0.070,0.075,0.080,0.085,0.090,0.095, \
      0.10, 0.11, 0.12, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,  \
      1,    1.1,  1.2,  1.5,  2,    2.5,  3,    3.5,  4,    4.5,  5,    5.5,  6,    6.5,  7,    7.5,  8,    8.5,  9,    9.5,   \
      10,   11,   12,   15,   20,   25,   30,   35,   40,   45,   50,   55,   60,   65,   70,   75,   80,   85,   90,   95,    \
      100,  110,  120,  150,  200,  250,  300,  350,  400,  450,  500,  550,  600,  650,  700,  750,  800,  850,  900,  950,   \
      1000, 1100, 1200, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 8500, 9000, 9500   \
    )
  }
  }
} else {
  unset logscale y
}

set grid xtics  lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 2 lc rgb '#ffddaa'
set grid mxtics
set grid ytics  lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 2 lc rgb '#ffddaa'

set xzeroaxis lt 1 lw 3 lc rgb '#ffddaa'

if (nfiles > 1) {
  set key left top
} else {
  unset key
}

load "${gplFile}"
quit
EOF

if [[ -s ${fullplotFile} ]]; then 
  convert ${fullplotFile} -resize '50%' ${plotFile}
  if [[ "/${show}" == "/YES" ]]; then display ${plotFile}; fi
  cat ${plotFile}
fi

rm -fv ${tmp}{-*,}.*
    
        
