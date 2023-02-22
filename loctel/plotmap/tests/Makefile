# Last edited on 2023-02-22 08:01:21 by stolfi

PROG := plotmap
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TOWNS := \
  bahia \
  campinas \
  itajobi
  
VERSIONS := \
  sq \
  hc

.PHONY:: run-all run-single  clean-single

all: clean run-all

clean::
	-/bin/rm -f out/*.eps
        
run-all: 
	for town in ${TOWNS} ; do \
          for vers in ${VERSIONS} ; do \
            ${MAKE} TOWN=$$town VERSION=$$vers clean-single run-single ; \
          done \
        done

# ######################################################################
# Begin section ${TOWN} ${VERSION}

ifneq "/${TOWN}" "/"
ifneq "/${VERSION}" "/"

PSVIEW := evince

ifeq "/${VERSION}" "/sq"
  CLIPOP := -showRect
else ifeq "/${VERSION}" "/hc"
  CLIPOP := -clip
endif

ifeq "/${TOWN}" "/bahia"
  RCIRCLE := 290
  RBALL := 410
  RECT := 6500 8100 6200 7800
else ifeq "/${TOWN}" "/campinas" 
  RCIRCLE := 212
  RBALL := 300
  RECT := 100 2300 400 2600
else ifeq "/${TOWN}" "/itajobi"
  RCIRCLE := 212
  RBALL := 300
  RECT := 100 2300 100 3100
endif

MAPFILE := in/${TOWN}-map.rnt
TUPSFILE := in/${TOWN}-${VERSION}.tups

TUPSDRAWFILE := out/${TOWN}-${VERSION}.opts
EPSFILE := out/${TOWN}-${VERSION}.eps
 
clean-single: 
	-/bin/rm -f ${EPSFILE} ${TUPSDRAWFILE}
 
run-single: ${EPSFILE}

${EPSFILE}: ${MAPFILE} ${TUPSDRAWFILE} ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
          -size 160 160 \
          -grid 1000 \
          ${CLIPOP} ${RECT} \
          ` cat ${TUPSDRAWFILE} ` \
          < ${MAPFILE} \
          > ${EPSFILE}
	${PSVIEW} ${EPSFILE}

${TUPSDRAWFILE}: ${TUPSFILE}
	cat ${TUPSFILE} \
          | gawk ' \
              /^ *([#]|$$)/{ next; } \
              /./{ \
                printf "-showCircle %d %d ${RCIRCLE} ", $$1, $$2; \
                printf "-showBall %d %d ${RBALL}\n", $$1, $$2; \
              } ' \
          > ${TUPSDRAWFILE}
endif
endif
# End section ${TOWN} ${VERSION}
# ######################################################################
