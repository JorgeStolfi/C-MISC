# Last edited on 2024-11-12 17:17:49 by stolfi

PROG := vigenere_wiki

all: ${PROG}
	${PROG}
        
${PROG}: ${PROG}.c
	gcc -c ${PROG}.c -o ${PROG}.o
	gcc -o ${PROG} ${PROG}.o -lm
