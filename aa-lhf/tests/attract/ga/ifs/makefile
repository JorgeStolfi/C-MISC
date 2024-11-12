#Last edited on 2007-04-16 21:16:04 by stolfi

CC=/usr/bin/gcc
IFLAGS := -I${STOLFIHOME}/lib/gp
CFLAGS := -g -O2 -Wall -ffloat-store -frounding-math
LIBS= -L${STOLFIHOME}/lib/gp -lgp -L/usr/X11R6/lib -lX11 -lm

all:	ia aa

aa:	ga.o aa.o
	${CC} ${CFLAGS} ${IFLAGS} -o aa ga.o aa.o ${LIBS}

ia:	ga.o ia.o
	${CC} ${CFLAGS} ${IFLAGS} -o ia ga.o ia.o ${LIBS}

ifs:	ifs.o
	${CC} ${CFLAGS} ${IFLAGS} -o ifs ifs.o ${LIBS}

clean:
	rm -f *.o ia aa ifs
