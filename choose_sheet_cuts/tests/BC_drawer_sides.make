# Last edited on 2020-03-08 11:50:26 by jstolfi
# Last edited on 2019-05-28 16:28:36 by jstolfi

PROG := choose_sheet_cuts
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: do-test

TESTNAME := BC_drawer_sides
 
OUTPREFIX := out/${TESTNAME}
INPREFIX := in/${TESTNAME}

OUTFILE := out/${TESTNAME}_layout.txt

do-test: do-clean ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
	  -stockSheets A 1 plywd 11 700 700 \
	  -stockSheets B 1 plywd 11 700 700 \
	  -stockSheets C 1 plywd 11 700 700 \
	  -stockSheets D 1 plywd 11 700 700 \
	  -stockSheets E 3 plywd 11 700 700 \
	  -stockSheets F 1 plywd 11 700 700 \
	  -stockSheets G 3 plywd 11 700 700 \
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
