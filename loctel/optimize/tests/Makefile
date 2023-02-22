# Last edited on 2023-02-22 09:32:26 by stolfi

PROG := optimize
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TOWNS := \
  bahia \
  campinas \
  itajobi

.PHONY:: run-all run-single clean-single

all: clean run-all

clean::
	-/bin/rm -f out/*.eps out/*.phn out/*.plt
        
run-all: 
	for town in ${TOWNS} ; do \
          ${MAKE} TOWN=$$town run-single clean-single ; \
        done

# ######################################################################
# Begin section ${TOWN}

ifneq "/${TOWN}" "/"

PSVIEW := evince
SHOWPGM := display

MAPNAME := in/${TOWN}-map
MAPFILE := ${MAPNAME}.rnt

USRNAME := in/${TOWN}-usr
USRFILE := ${URNAME}.vdt

OUTPREFIX := out/${TOWN}

# Phone coordinate files:
PHNFILES := \
  ${OUTPREFIX}-strip.phn \
  ${OUTPREFIX}-hexag.phn \
  ${OUTPREFIX}-greed.phn
  
EPSFILES := ${subst .phn,.eps,${PHNFILES}}

PLTOPSFILES := ${subst .plt,.eps,${PHNFILES}}

clean-single: rm -f ${EPSFILES} ${PHNFILES} ${PLTOPSFILES}
  
run-single: ${EPSFILES}

# ANATEL's maximum walking distance:
MAXDIST := 400

# ${MAXDIST} plus about half the maximum edge length:
MAXDISTPLUS := 410

# Distance at which capture factor falls to 0.5:
PARTDIST := 200

# Honeycomb radius: ${MAXDIST} times sqrt(0.5) minus quantization error:
HEXRADIUS := 290

ifeq "${TOWN}" "bahia" 
  STRIPWIDTHS := 40 0.25 1.50
else
  STRIPWIDTHS := 20 0.25 2.00
endif

${PHNFILES}: ${MAPFILE} in/${TOWN}-usr.vdt ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
          -mapName ${MAPNAME} \
          -demName ${USRNAME} \
          -maxDist ${MAXDIST} \
          -partDist ${PARTDIST} \
          -stripWidths ${STRIPWIDTHS} \
          -outName ${OUTPREFIX}
	wc -l ${PHNFILES}
        
PLOTPROG := ../../plotmap/plotmap

${OUTPREFIX}-%.plt: ${OUTPREFIX}-%.phn
	cat ${OUTPREFIX}-$*.phn \
          | gawk \
              -v hr=${HEXRADIUS} \
              -v md=${MAXDISTPLUS} \
              ' /phones =/ { next; } \
                //{ \
                    printf "-showCircle %d %d %d  -showBall %d %d %d\n", \
                      $$3, $$4, hr, $$3, $$4, md; \
                  } ' \
          > ${OUTPREFIX}-$*.plt

${OUTPREFIX}-%.eps: ${MAPFILE} ${OUTPREFIX}-%.plt ${PLOTPROG}
	${PLOTPROG} \
          -size 160 160 \
          -grid 1000 \
          ` cat ${OUTPREFIX}-$*.plt ` \
          < ${MAPFILE} > ${OUTPREFIX}-$*.eps
	${PSVIEW} ${OUTPREFIX}-$*.eps

endif
# End section ${TOWN}
# ######################################################################
