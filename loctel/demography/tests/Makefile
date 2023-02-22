# Last edited on 2023-02-22 08:00:02 by stolfi

PROG := demography
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TOWNS := \
  bahia \
  campinas \
  itajobi

.PHONY:: run-all run-single clean-single 

all: clean run-all

clean:: 
	-/bin/rm -f out/*.pgm out/*.vdt
        
run-all: 
	for town in ${TOWNS} ; do \
          ${MAKE} TOWN=$$town clean-single run-single ; \
        done

# ######################################################################
# Begin section ${TOWN}

ifneq "/${TOWN}" "/"

PSVIEW := evince
SHOWPGM := display

OUTPREFIX := out/${TOWN}

PGMFILES := \
  ${OUTPREFIX}-std.pgm \
  ${OUTPREFIX}-pop.pgm

VDTFILES := \
  ${OUTPREFIX}-usr.vdt \
  ${OUTPREFIX}-vms.vdt

clean-single: 
	-/bin/rm -f ${OUTPREFIX}-*.pgm ${OUTPREFIX}-*.vdt 
  
run-single: ${PGMFILES}

${PGMFILES}: in/${TOWN}-map.rnt ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
          -mapName in/${TOWN}-map \
          -step 25 -radius 750 -totPop 50000 \
          -outName ${OUTPREFIX}
	${SHOWPGM} ${PGMFILES}          

endif
# End section ${TOWN}
# ######################################################################
