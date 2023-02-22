#! /bin/bash -f
# Last edited on 2023-02-12 10:05:16 by stolfi

eps_files=( "$@" )
pdf_files=()

for efile in ${eps_files[@]} ; do
  pfile="${efile%.*}.pdf"
  echo "converting ${efile} to ${pfile}" 1>&2
  echo "quit" \
    | gs -sOutputFile="-" -sDEVICE=pdfwrite -g6253x4576 -r72 -sPAPERSIZE=a4 -dEPSFitPage ${efile} \
    > ${pfile} 
  pdf_files+=( ${pfile} )
done 
for ff in ${pdf_files[@]} ; do 
  { evince ${ff} & }
done
display Q.png

# -dFIXEDMEDIA
# -dFIXEDMEDIA  -sDEVICE=pdfwrite -sOutputFile=\temp\out.pdf -c "<</PageSize [968 1184]>> setpagedevice -20 -50 translate" -f d:\temp\ori.eps

# gs -o myfile.pdf -sDEVICE=pdfwrite -g5775x6207 -dPDFFitPage myfile.ps
# -g... gives the medium size in pixels.
# An A4 page has a dimension of 595x842pt (PostScript points).
# 1 Inch is the same as 72 PostScript points.
# Ghostscript internally by default computes with a resolution of 720 pixels per inch when it comes to PDF output.
# Hence for PDF output 595x842pt == 5950x8420px.
# Hence for your case in question 8.02x8.62Inches ≈≈ 5775x6207px.
