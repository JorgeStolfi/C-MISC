# Last edited on 2020-03-08 11:45:03 by jstolfi
# Last edited on 2019-05-28 16:28:36 by jstolfi

PROG := choose_sheet_cuts
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: do-test

TESTNAME := BC_drawers
 
OUTPREFIX := out/${TESTNAME}
INPREFIX := in/${TESTNAME}

OUTFILE := out/${TESTNAME}_layout.txt

do-test: do-clean ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
	  -stockSheets D 1 plywd 11 1050  550 \
	  -stockSheets B 1 plywd 11 1250  350 \
	  -stockSheets C 1 plywd 11  700  550 \
	  -stockSheets F 1 plywd 11 1550  300 \
	  -stockSheets A 3 plywd 11  450  600 \
	  -stockSheets E 1 plywd 11 2200 1000 \
	  -stockSheets X 3 plywd 11 2200 1600 \
	  -scrapMargin 1 \
	  -cutWidth 4 \
	  -outPrefix ${OUTPREFIX} \
          -inPrefix ${INPREFIX}
	for ff in ${OUTPREFIX}_*.eps ; do \
	  pfile="$${ff%.*}.pdf"; \
          echo "converting $${ff} to $${pfile}" 1>&2 ; \
	  ps2pdf -sPAPERSIZE=letter -dFIXEDMEDIA -dEPSFitPage $${ff} $${pfile} ; \
	done
	atril ${OUTPREFIX}_*.eps
	atril ${OUTPREFIX}_*.pdf

do-clean:
	rm -f ${OUTPREFIX}_*.eps ${OUTPREFIX}_*.pdf
