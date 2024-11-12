CC=/usr/bin/gcc
CFLAGS := -g -O2 -Wall -ffloat-store -frounding-math
LIBS= -L/n/X11R6/lib -lX11 -lm

all:	ia aa

aa:	ga.o aa.o f.o
	${CC} ${CFLAGS} ${IFLAGS} -o aa ga.o aa.o f.o ${LIBS}

ia:	ga.o ia.o f.o
	${CC} ${CFLAGS} ${IFLAGS} -o ia ga.o ia.o f.o ${LIBS}

henon:	henon.o
	${CC} ${CFLAGS} ${IFLAGS} -o henon henon.o ${LIBS}
clean:
	rm -f *.o ia aa
