#! /bin/bash

TAG="$1"; shift
NTERMS="$1"; shift
NOMIDDLE="$1"; shift
RANGES=(`echo "$1" | sed -e 's:[/]: :g'` ); shift

PROGDIR=".."
PROG="gauss_tripart_opt"
TMP="/tmp/$$"
rfile="${TMP}.args" # File with the range arguments. 


for k in `count 0 ${NTERMS}` ; do
  rark=(`echo "${RANGES[k]}" | sed -e 's:[,]: :g'` )
  if [[ $k -ne ${rark[0]} ]]; then
    echo "** term index mismatch" 1>&2; exit 1
  fi
  if [[ $k -gt 0 ]]; then
    printf "%s " "-rangeA ${k} ${rark[1]} ${rark[2]}" >> ${rfile}
  fi
  printf "%s\n" "-rangeD ${k} ${rark[3]} ${rark[4]} " >> ${rfile}
done

${PROGDIR}/${PROG} \
  -outPrefix out/funcs_${TAG} \
  -nTerms ${NTERMS} \
  -noMiddle ${NOMIDDLE} \
  -nSamples 1001 \
  -minD 0.200 \
  -maxA 3.000 \
  -minSepA 0.25 \
  `cat ${rfile}`
  
rm -f ${TMP}-*

