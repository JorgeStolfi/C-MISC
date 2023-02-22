# Last edited on 2023-02-09 08:09:20 by stolfi

PROG := bmix
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

IN_FILES := data/afile.txt data/bfile.txt
KEY_FILE := data/xfile.txt
OUT_FILES := out/amfile.txt out/bmfile.txt
OUT2_FILES := out/aufile.txt out/bufile.txt

all: out/amfile.txt out/aufile.txt

# Mixing two files:
out/amfile.txt: ${PROGDIR}/${PROG} ${IN_FILES}
	${PROGDIR}/${PROG} \
            -verbose \
            ${IN_FILES} ${KEY_FILE} ${OUT_FILES}
	ls -l ${IN_FILES} ${KEY_FILE} ${OUT_FILES}
	strings out/amfile.txt

# Un-mixing the two files:
out/aufile.txt: ${PROGDIR}/${PROG} ${OUT_FILES} ${KEY_FILE}
	${PROGDIR}/${PROG} \
            -verbose \
            ${OUT_FILES} ${KEY_FILE} ${OUT2_FILES} 
	ls -l ${OUT_FILES} ${KEY_FILE} ${OUT2_FILES}
	prdiff data/afile.txt out/aufile.txt
	prdiff data/bfile.txt out/bufile.txt
