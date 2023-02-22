# Last edited on 2023-02-08 22:44:53 by stolfi

PROG := bxor
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: out/abfile.txt out/abbfile.txt out/aafile.txt

out:
	mkdir -p out

out/abfile.txt: ${PROGDIR}/${PROG} data/afile.txt data/bfile.txt out 
	${PROGDIR}/${PROG} \
            -verbose \
            data/afile.txt data/bfile.txt \
          > out/abfile.txt
	ls -l data/afile.txt data/bfile.txt out/abfile.txt
	strings out/abfile.txt

out/abbfile.txt: ${PROGDIR}/${PROG} out/abfile.txt out data/bfile.txt
	${PROGDIR}/${PROG} \
            -verbose \
            out/abfile.txt data/bfile.txt \
          > out/abbfile.txt
	ls -l out/abfile.txt data/bfile.txt out/abbfile.txt
	prdiff data/afile.txt out/abbfile.txt

out/aafile.txt: ${PROGDIR}/${PROG} data/afile.txt out
	${PROGDIR}/${PROG} \
            -verbose \
            data/afile.txt data/afile.txt \
          > out/aafile.txt
	ls -l data/afile.txt data/afile.txt out/aafile.txt
	od out/aafile.txt
