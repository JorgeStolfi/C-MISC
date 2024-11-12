# -*- makefile -*-
# $Id: Makefile,v 1.7 1997/06/19 19:09:58 pmartin Exp $

BIM_FLAGS := -DNO_ERROR_HANDLING -DSTATS -DCONCISE
CXXFLAGS += ${BIM_FLAGS}

BINARIES = bim_hi bim_lo bim_fifo bim_lifo

MAIN = bim_run.cc

# pick options for the binaries
bim_lo_opts=-DBUCKET -DLO_SELECT 
bim_hi_opts=-DBUCKET -DHI_SELECT 
bim_fifo_opts= -DNO_BUCKET -DQUEUE -DFIFO 
bim_fifo_b_opts=-DBUCKET -DQUEUE -DFIFO 
bim_lifo_opts=-DNO_BUCKET -DQUEUE -DLIFO 
bim_lifo_b_opts=-DBUCKET -DQUEUE -DLIFO 

include ../../make.config

