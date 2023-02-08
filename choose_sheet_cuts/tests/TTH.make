# Last edited on 2022-04-20 00:17:33 by stolfi
# Last edited on 2019-05-28 16:28:36 by jstolfi

PROG := choose_sheet_cuts
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: do-test

TESTNAME := TTH
 
OUTPREFIX := out/${TESTNAME}
INPREFIX := in/${TESTNAME}

OUTFILE := out/${TESTNAME}_layout.txt

INPLATES := ${INPREFIX}_plates.txt
INGROUP := ${INPREFIX}_grouping.txt

SHOW_EPS := NO
SHOW_PDF := YES

do-test: do-clean ${PROGDIR}/${PROG} ${INPLATES} ${INFGROUP}
	${PROGDIR}/${PROG} \
	  -stockSheets A 1 ac.m    3.0 1000  500 \
	  -stockSheets B 1 psai.o  3.0 1000  500 \
	  -stockSheets C 1 psai.o  5.0 1000  500 \
	  -stockSheets D 1 ps.d    1.7  600  300 \
	  -stockSheets E 1 ps.m    2.0 1000  495 \
           \
	  -scrapMargin 2 \
	  -cutWidth 3 \
	  -outPrefix ${OUTPREFIX} \
          -inPrefix ${INPREFIX}
	ls -l ${OUTPREFIX}_*.eps
	if [[ "/${SHOW_EPS}" == "/YES" ]]; then \
	  for ff in ${OUTPREFIX}_*.eps ; do \
            ( evince $${ff} & ) ; \
          done ; \
        fi
	if [[ "/${SHOW_PDF}" == "/YES" ]]; then \
          convert_eps_to_pdf.sh ${OUTPREFIX}_*.eps ; \
        fi

do-clean:
	rm -f ${OUTPREFIX}_*.eps ${OUTPREFIX}_*.pdf
