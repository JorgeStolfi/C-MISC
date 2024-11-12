# Last edited on 2024-11-11 23:35:40 by stolfi

PROG := avgdet

JS_LIBS = \
  libgeo.a \
  libjs.a
  
all: build runit
  
include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make

runit: ${PROG}
	${PROG}
