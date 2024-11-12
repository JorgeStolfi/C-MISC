# -*- makefile -*-

C_LIBS += -lm
OBJS = random.o
BINARIES = rbg1 rbg2

include ../../make.config

rbg1: rbg1.o
rbg2: rbg2.o

ifndef ndep
-include .rbg1.d .rbg2.d
endif
