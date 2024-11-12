CC=gcc
CFLAGS= -g -O2 -Wall -I$(HOME)/lib/gp
LIBS= -L$(HOME)/lib/gp -lgp -L/usr/X11R6/lib -lX11 -lm

all:	ia aa

aa:	ga.o aa.o f.o
	$(CC) $(CFLAGS) -o aa ga.o aa.o f.o $(LIBS)

ia:	ga.o ia.o f.o
	$(CC) $(CFLAGS) -o ia ga.o ia.o f.o $(LIBS)

henon:	henon.o
	$(CC) $(CFLAGS) -o henon henon.o $(LIBS)
clean:
	rm -f *.o ia aa
