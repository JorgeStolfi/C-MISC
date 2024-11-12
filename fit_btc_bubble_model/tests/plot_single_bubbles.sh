#! /bin/bash
# Last edited on 2015-05-04 14:16:39 by stolfilocal

# Plots individual bubbles, lin and log scale

rundate="$1"; shift; # Nominal date of experiment.
hrad="$1"; shift;    # Radius of price smoothing window.
nb="$1"; shift;      # Number of bubbles in parameter file.
inFile="$1"; shift;  # Filename with given price, model price, bubble components.

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

# Get the bubbles colors {bcols}, tags {btags} and number of bubbles {nb}: 
btags=( `cat ${parmsFile} | gawk '/^ *[-+]?[.0-9]+[ ]+20[01][0-9]-/{ print $8; }'` )
bcols=( `cat ${parmsFile} | gawk '/^ *[-+]?[.0-9]+[ ]+20[01][0-9]-/{ print $9; }'` )
if [[ ( 10#${nb} -ne ${#btags[@]} ) || (10#${nb} -ne ${#bcols[@]} ) ]]; then
  echo "** wrong number of bubbles nb = ${nb} found ${#btags[@]} ${#bcols[@]}" 1>&2
  exit 1
fi

tmp="/tmp/$$"

ib=0
while [[ ${ib} -lt 10#${nb} ]]; do

  dmin="2010-06-20"    # Min date for plot.
  dmax="2015-05-10"    # Max date for plot.

  make_jpg="YES"
  show_jpg="NO"
  show_png="NO"

  ibx=`printf '%02d' "${ib}"`
  
  for logScale in 0 1 ; do 

    # Extract the bubble, and the relevant part of the price and residual
    tmpFile="${tmp}.dat"
    gawk \
        -v ib=${ib} \
        ' /^ *([\#]|$)/ { print; next; } 
          /^ *20[01][0-9][-]/ { 
            dt = $1; tm = $2;
            ap = $4; # Actual price.
            mp = $6; # Modeled price.
            rs = $8; # Residual {ap-mp}.
            kb = 12 + 2*ib; # Index of bubble price field.
            bp = $(kb); # Bubble price.
            
            printf "%s %s ", dt, tm;
            if (bp + 0 == 0) { 
              pap = 0.0; # Actual price to plot.
              pmp = 0.0; # Model price to plot.
              prs = 0.0; # Residual price to plot.
              prb = 0.0; # Residual plus bubble.
            } else {
              pap = ap;  # Actual price to plot.
              pmp = mp;  # Model price to plot.
              prs = rs;  # Residual price to plot.
              prb = rs + bp; # Residual plus bubble.
            }
            printf " %18.5f %18.5f %18.5f %18.5f %18.5f\n", pap, pmp, prs, bp, prb;
            next;
          }
          // { printf "** bad format\n" > "/dev/stderr"; exit(1); }
        ' \
      < ${inFile} \
      > ${tmpFile}

    if [[ ${logScale} -gt 0 ]]; then
      # Plot the actual and model price too, not the residuals:
      prArgs=( "${tmpFile}" 3 "Price" 1.000 778899 "${tmpFile}" 4 "Model" 1.000 228822 )
      vmin="0.04"
      vmax="1600.0"
    else
      # Plot the bubble and the residual that the bubble is supposed to model:
      # prArgs=( "${tmpFile}" 7 "Target" 1.000 ee2200 )
      
      # Plot only the bubble:
      prArgs=(  )
      vmin="0.00"
      vmax="0.00"
    fi
    
    plotPref="${outPref}-bub${ibx}-log${logScale}"; # Plot file prefix.

    pngFile="${plotPref}.png"
    jpgFile="${plotPref}.jpg"
    rm -f ${pngFile} ${jpgFile}
    
    
    plot_prices.sh \
      "Price bubble ${ibx} (${btags[${ib}]}) - ${dmin} to ${dmax}${wintit}" \
      180 6 \
      "${dmin}" "${dmax}" \
      "${vmin}" "${vmax}" \
      ${logScale} \
      "${show_png}" \
      ${prArgs[@]} \
      "${tmpFile}" 6 "${btags[${ib}]}" 1.000 "${bcols[${ib}]}" \
      > ${pngFile}
      
    rm -fv ${tmpFile}

    if [[ ( -s ${pngFile} ) && ( "/${make_jpg}" == "/YES") ]]; then
      convert ${pngFile} -quality 85 -resize '600x' ${jpgFile}
      ls -l ${pngFile} ${jpgFile}
      if [[ "/${show_jpg}" == "/YES" ]]; then display ${jpgFile}; fi
    fi
  done
  
  ib=$(( ${ib} + 1 ))
done
  
# rm -fv ${tmp}-rem-*.txt
