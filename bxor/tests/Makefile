# Last edited on 2013-10-31 02:27:22 by stolfilocal

PROG := bxor
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: out/abfile.txt out/abbfile.txt out/aafile.txt

out:
	mkdir -p out

out/abfile.txt: ${PROGDIR}/${PROG} afile.txt bfile.txt out 
	${PROGDIR}/${PROG} \
            -verbose \
            afile.txt bfile.txt \
          > out/abfile.txt
	ls -l afile.txt bfile.txt out/abfile.txt
	strings out/abfile.txt

out/abbfile.txt: ${PROGDIR}/${PROG} afile.txt out/abfile.txt out
	${PROGDIR}/${PROG} \
            -verbose \
            out/abfile.txt bfile.txt \
          > out/abbfile.txt
	ls -l out/abfile.txt bfile.txt out/abbfile.txt
	prdiff afile.txt out/abbfile.txt

out/aafile.txt: ${PROGDIR}/${PROG} afile.txt out
	${PROGDIR}/${PROG} \
            -verbose \
            afile.txt afile.txt \
          > out/aafile.txt
	ls -l afile.txt afile.txt out/aafile.txt
	od out/aafile.txt
