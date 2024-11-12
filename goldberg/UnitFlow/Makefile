# -*- makefile -*-

# This makefile recursively processes the subdirectories listed in the 
# variable SUBDIRS. To add a new generator, drop it in a subdirectory 
# (with its makefile) and add the subdirectory to SUBDIRS here.
# The targets all, install, and clean are expected in the Makefile in 
# the subdirectory.
SUBDIRS=pr bfs dfs

all:

%: 
	for i in ${SUBDIRS}; do cd $$i; ${MAKE} $@; cd ..; done
