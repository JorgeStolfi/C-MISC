# Last edited on 2023-02-03 23:22:34 by stolfi

PROGS := \
  SOPlotFunction \
  SOMakeProcFunction \
  SOMakeWaveBasis

JS_LIBS := \
  libdygrid.a \
  libeps.a \
  libgeo.a \
  libjs.a

LIB_HFILES_NOT := \
  SOBezFunction.h \
  SOSpline.h \
  SOFuncMap.h \
  SOFunctionOrtho.h

LIB_CFILES_NOT := \
  SOBezFunction.c \
  SOFuncMap.c \
  SOSpline.c \
  SOFunctionOrtho.c
  
PROGS_LATER := \
  SOComputeApprox \
  SOTestTentFunction \
  SOTestIntegral \
  SOComputeErrorMap \
  SOUniSolve \
  SOTestGrid \
  SOMakeGrid \
  SOFindTentBasis \
  SOComputeMatrices
 
PROGS_NOT_YET := \
  SOCompareMatrices \
  SOComputeMatrices \
  SOExtractBasis \
  SOMakePolyBasis

IGNORE := \
  ${LIB_HFILES_NOT} \
  ${LIB_CFILES_NOT} \
  ${addsuffix .c,${PROGS_LATER}} \
  ${addsuffix .c,${PROGS_NOT_YET}}

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
