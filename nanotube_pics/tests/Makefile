# Last edited on 2023-02-03 19:46:53 by stolfi

PROGDIR := ..
PROG := nanotube_pics

all: ${PROGDIR}/${PROG}
	rm -f out/*.eps out/*.pdf
	${PROGDIR}/${PROG}
	for psfile in out/*.eps ; do \
	  pdffile="$${psfile%%.eps}.pdf" ; \
	  echo "$${psfile} -> $${pdffile}" ; \
	  ps2pdf \
	    -dEPSCrop \
	    -dPDFSETTINGS=/prepress \
	    -dEmbedAllFonts=true \
	     $${psfile} $${pdffile} ; \
	  pdffonts $${pdffile} ; \
	done
	atril out/*.pdf
