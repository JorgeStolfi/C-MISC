# Last edited on 2012-12-15 12:03:32 by stolfilocal

PROG := raut_build
PROGDIR := ../..

all: run

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

run: ${PROGDIR}/${PROG}
	cat in/portu-ini.dic \
        | ( export LC_ALL=C; sort ) \
        | ${PROGDIR}/${PROG} \
              -crunch 10000 \
              -plot   out/portu-ini.plt 2000 \
              -add    -
	../plot_growth.sh out/portu-ini
