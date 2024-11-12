# -*- makefile -*-
NO_INSTALL = 1
BINARIES = pow
MAIN = pow.c
C_LIBS = -lm
install: ${BINARIES}
clean:
	rm -f ${BINARIES}
realclean: clean
	rm -f $(BINARIES:%=.%.d)
uninstall: 
	rm -f ${BINARIES}

include ../make.config
