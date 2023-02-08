# Last edited on 2010-05-21 12:42:07 by stolfi

PROG := nclassif

JS_LIBS := \
  libimg.a \
  libjs.a

all: build install

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make

