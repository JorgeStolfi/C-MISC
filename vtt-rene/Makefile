# Makefile by J. Stolfi
# Last edited on 2008-01-14 22:41:51 by stolfi

PROG := vtt

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make

CURVERSION := 1.7e1

OLDVERSION := 1.7

compare: ${PROG}.diffs

${PROG}.diffs: ${PROG}.c
	-diff -Bb -C2 ${OLDVERSION}/${PROG}.c ${PROG}.c > ${PROG}-${OLDVERSION}-${CURVERSION}.diffs

# ----------------------------------------------------------------------
# "make export-archive" is an alternative to the standard "make archive":

archive: ${SOURCE} ${PROG}.diffs
	tar -cvzf vtt-${CURVERSION}.tgz \
          README.DOC ${PROG}.doc ${PROG}.c Makefile \
          ${PROG}.diffs \
          tests/Makefile tests/do-all-tests \
          tests/text*.evt tests/cmp*.evt
