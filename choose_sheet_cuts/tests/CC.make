# Last edited on 2020-11-06 22:53:22 by jstolfi
# Last edited on 2019-05-28 16:28:36 by jstolfi

PROG := choose_sheet_cuts
PROGDIR = ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: do-test

TESTNAME := CC
 
OUTPREFIX := out/${TESTNAME}
INPREFIX := in/${TESTNAME}

OUTFILE := out/${TESTNAME}_layout.txt

SHOW_EPS := NO
SHOW_PDF := YES

do-test: do-clean ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
	  -stockSheets H 1 plywd  4.5 2200 1600 \
	  -stockSheets D 1 plywd 11.0  570 1050 \
	  -stockSheets F 1 plywd 11.0  400  880 \
	  -stockSheets J 1 plywd 11.0 2200 1600 \
 	  -stockSheets E 1 plywd 11.0 1630  360 \
	  -stockSheets G 1 plywd 11.0  370  560 \
          \
	  -stockSheets X 1 plywd  4.5 2200 1600 \
	  -stockSheets Y 3 plywd 11.0 2200 1600 \
	  -stockSheets Z 3 plywd 15.5 2200 1600 \
           \
	  -scrapMargin 5 \
	  -cutWidth 5 \
	  -outPrefix ${OUTPREFIX} \
          -inPrefix ${INPREFIX}
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
