# Last edited on 2013-10-31 20:19:56 by stolfilocal

PROG := bmix
PROGDIR := ../..

include ${PROGDIR}/GENERIC-PROGS-TEST.make

all: out/amfile.txt out/aufile.txt

out:
	mkdir -p out

# Mixing two files:
out/amfile.txt: ${PROGDIR}/${PROG} afile.txt bfile.txt xfile.txt out 
	${PROGDIR}/${PROG} \
            -verbose \
            afile.txt bfile.txt \
            xfile.txt \
            out/amfile.txt out/bmfile.txt
	ls -l afile.txt bfile.txt xfile.txt out/amfile.txt out/bmfile.txt
	strings out/amfile.txt
	strings out/bmfile.txt
	diff out/amfile.txt out-check/amfile.txt
	diff out/bmfile.txt out-check/bmfile.txt

# Un-mixing the two files:
out/aufile.txt: ${PROGDIR}/${PROG} out/amfile.txt out/bmfile.txt xfile.txt out
	${PROGDIR}/${PROG} \
            -verbose \
            out/amfile.txt out/bmfile.txt \
            xfile.txt \
            out/aufile.txt out/bufile.txt
	ls -l out/{a,b}mfile.txt xfile.txt out/{a,b}ufile.txt
	diff afile.txt out/aufile.txt
	diff bfile.txt out/bufile.txt
