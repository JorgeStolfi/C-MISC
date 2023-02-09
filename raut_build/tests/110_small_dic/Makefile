# Last edited on 2012-12-15 12:03:46 by stolfilocal

PROG := raut_build
PROGDIR := ../..

all: run

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

run: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -plot   out/small-ini.plt 5 \
          -add    in/small-ini.dic \
          -remove in/small-sub.dic \
          -add    in/small-add.dic \
          -modify in/small-mod.dic \
          -check  in/small-chk.dic out/small-bad.dic \
          -list - \
          -dump   out/small.dmp
	../plot_growth.sh out/small-ini
