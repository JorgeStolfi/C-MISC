#! /bin/bash
# Last edited on 2015-05-02 23:05:16 by stolfilocal

# Plots all bubbles together

rundate="$1"; shift; # Nominal date of experiment.
hrad="$1"; shift;    # Radius of price smoothing window.
nb="$1"; shift;      # Number of bubbles in parameter file.
inFile="$1"; shift;  # Filename with given price, model price, component bubbles.

dmin="2010-06-20"    # Min date for plot.
dmax="2015-05-10"    # Max date for plot.

make_jpg="YES"
show_jpg="NO"

trap "exit 1" SIGINT SIGTERM

dataDir="data"
outDir="out"

parmsFile="${dataDir}/${rundate}-bubble-${nb}.parms"

outPref="${outDir}/${rundate}-bubble-${nb}-fit-sm${hrad}"

# Filtering note on title:
if [[ 10#${hrad} -eq 0 ]]; then
  wintit=""
else
  hwid=$(( 2 * 10#${hrad} + 1 )); # Total width of smoothing window.
  wintit=" (smoothed with ${hwid}-day Hann window)"
fi

logScale=1

plotPref="${outPref}"; # Plot file prefix.
pngFile="${plotPref}.png"
jpgFile="${plotPref}.jpg"
rm -f ${pngFile} ${jpgFile}

plot_prices.sh \
  "Price bubbles - ${dmin} to ${dmax}${wintit}" \
  180 6 \
  "${dmin}" "${dmax}" \
  "0.04" "1600.0" \
  ${logScale} \
  "YES" \
  "${inFile}"  4 "Price"    1.000 778899 \
  "${inFile}"  6 "Model"    1.000 228822 \
  "${inFile}" 10 "Ratio"    1.000 882222 \
  ` cat ${parmsFile} \
      | gawk \
          -v hrad=${hrad} \
          -v ifi=${inFile} \
          ' BEGIN { kf = 12; } \
            /^ *[-+]?[.0-9]+[ ]+20[01][0-9]-[01][0-9]-[0-3][0-9]/ {
              printf "%s %d %s 1.000 %s\n", ifi, kf, $8, $9; 
              kf += 2;
            } \
          ' \
  ` \
  > ${pngFile}

# To plot the residual:
# "${inFile}"  8 "Residual" 1.000 ff9900 \
#   

if [[ ( -s ${pngFile} ) && ( "/${make_jpg}" == "/YES") ]]; then
  convert ${pngFile} -quality 85 -resize '600x' ${jpgFile}
  ls -l ${pngFile} ${jpgFile}
  if [[ "/${show_jpg}" == "/YES" ]]; then display ${jpgFile}; fi
fi

  
# rm -fv ${tmp}-rem-*.txt
