# -*- makefile -*-
# $Id: Makefile,v 1.7 1997/08/15 17:51:08 pmartin Exp $

UF_FLAGS := -DSTATS -DNO_ERROR_HANDLING -DPHASE2 -DCONCISE -DEXCESSES
CXXFLAGS += ${UF_FLAGS}

BINARIES = uf_lo uf_hi uf_fifo uf_lifo

MAIN = uf_run.cc

# pick options for the binaries
uf_lo_opts=-DBUCKET -DLO_SELECT 
uf_hi_opts=-DBUCKET -DHI_SELECT 
uf_fifo_opts= -DNO_BUCKET -DQUEUE -DFIFO 
uf_fifo_b_opts=-DBUCKET -DQUEUE -DFIFO 
uf_lifo_opts=-DNO_BUCKET -DQUEUE -DLIFO 
uf_lifo_b_opts=-DBUCKET -DQUEUE -DLIFO 

include ../../make.config

