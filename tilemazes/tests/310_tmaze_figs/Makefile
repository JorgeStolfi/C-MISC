# Last edited on 2009-11-08 23:06:20 by stolfi

PROGDIR := ../..
PROG := tmaze_figs

.PHONY:: all-figs show-all-figs

all: all-figs show-all-figs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all-figs: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG}

show-all-figs: 
	display -title '%d/%f' `export LC_ALL=C; ls out/*/*.eps | sort`
        
