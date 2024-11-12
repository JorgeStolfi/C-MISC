CC=/usr/bin/gcc
IFLAGS := -I${STOLFIHOME}/lib/gp
CFLAGS := -g -O2 -Wall -ffloat-store -frounding-math
LIBS= -lgp -L/n/X11R6/lib -lX11 -lm

all:	ia aa

aa:	aa.c aa.h
	${CC} ${CFLAGS} ${IFLAGS} -c aa.c

ia:	ia.c 
	${CC} ${CFLAGS} ${IFLAGS} -c ia.c

clean:
	rm -f *.o
