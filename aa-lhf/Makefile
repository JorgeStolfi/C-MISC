# Last edited on 2023-01-14 12:17:21 by stolfi

# NONSTANDARD MAKEFILE

# BEWARE: The sources named "fn*.c" all implement the same three functions,
# so only one of those sources should be linked into the program. 

X := -L/usr/lib -lX11

# acc
# CC := acc
# STRICT := -Xc -vc
# DEBUG := 

# /usr/bin/gcc
STRICT := -Wall
DEBUG := -g

IFLAGS := -I. -I${STOLFIHOME}/lib/include
CFLAGS := ${STRICT} ${DEBUG} -ffloat-store -frounding-math
LIBS := $X -lm

CC := /usr/bin/gcc ${CFLAGS} ${IFLAGS} 

PROGS := \
  aer-aa \
  aer-ia \
  acheck-ia \
  acheck-aa \
  qei-aa \
  qei-ia 

FNOBJS := \
  fnmul.o \
  fnsqr.o \
  fntanh.o \
  fncoeur.o \
  fncharges.o
  
# GPO := nogp.o
GPO := dv.o gp.o
  
FN := fncharges

all: build install

build: ${PROGS} ${FNOBJS}

install: 

uninstall:

clean::
	rm -f *.o ${PROGS}

actions:
	${MAKE} ${ACTIONS} 

aer-aa: Makefile aer.o aalhf.o ${GPO} ${FN}.o
	${CC} ${CFLAGS} ${IFLAGS} -o aer-aa aer.o aalhf.o ${GPO} ${FN}.o ${LIBS}

aer-ia: Makefile aer.o ialhf.o ${GPO} ${FN}.o
	${CC} ${CFLAGS} ${IFLAGS} -o aer-ia aer.o ialhf.o ${GPO} ${FN}.o ${LIBS}


acheck-aa: Makefile acheck.o aalhf.o ${FN}.o
	${CC} ${CFLAGS} ${IFLAGS} -o acheck-aa acheck.o aalhf.o ${FN}.o ${LIBS}

acheck-ia: Makefile acheck.o ialhf.o ${FN}.o
	${CC} ${CFLAGS} ${IFLAGS} -o acheck-ia acheck.o ialhf.o ${FN}.o ${LIBS}
        
qei-aa: Makefile qei.o aalhf.o ${GPO} 
	${CC} ${CFLAGS} ${IFLAGS} -o qei-aa qei.o aalhf.o ${GPO} ${LIBS}
        
qei-ia: Makefile qei.o ialhf.o ${GPO} 
	${CC} ${CFLAGS} ${IFLAGS} -o qei-ia qei.o ialhf.o ${GPO} ${LIBS}
        
aer.o:: arith.h gp.h fn.h

acheck.o:: arith.h fn.h

qei.o:: arith.h gp.h 

aalhf.o:: arith.h

ialhf.o:: arith.h

nogp.o:: nogp.h

dv.o:: gp.h

gp.o:: gp.h

${FN}.o:: fn.h arith.h


