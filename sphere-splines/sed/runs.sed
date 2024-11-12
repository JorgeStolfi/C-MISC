s@^BIN = /home/staff/anamaria/${PLATFORM}/bin@PROGDIR := ../../progs@
/JSINC/i\
DATADIR := ../../data\

/^MAKEFILE *=/d
s@{BIN}@{PROGDIR}@g
s@[$]{MAKEFILE}@Makefile@g
s@ *-f Makefile@ @g
s@[.][.]/converted@\${PROGDIR}@g
s@PROG=@@
s@prog-single@@
