# Last edited on 2012-12-20 20:44:19 by stolfilocal

PROG := overlap_comb
PROGDIR := ..
        
include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: ${PROGDIR}/${PROG}
	-rm -f out/comb-*.txt
	${PROGDIR}/${PROG} 1 61 1 10
	for f in `ls out/comb-*.txt | sort -r` ; do \
          ( plot_comb_overlap.sh $$f & ); \
        done
	xmessage "OK"
