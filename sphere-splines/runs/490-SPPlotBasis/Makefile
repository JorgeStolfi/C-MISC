# Last edited on 2023-02-04 07:02:58 by stolfi

PROGDIR := ../../progs
DATADIR := ../../data

# BASIS

BASES := \
  pw/icosa/raw/nh-d6c1-NO

BASES_LATER := \
  gl/harm/d005 \
  gl/harm/d006 \
  gl/harm/d017 \
  pw/icosa/raw/h-d7c1-OLD \
  pw/icosa/raw/h-d7c1 \
  pw/icosa/raw/h-d7c1-NO \
  pw/icosa/raw/h-d7c1-OLD-NO \
  pw/icosa/raw/nh-d6c1-OLD-NO \
  pw/icosa/raw/nh-d6c1
  
BASFILES := ${addprefix ${DATADIR}/sp/,${addsuffix /b.bas,${BASES}}}
PSFILES  := ${addprefix ${DATADIR}/sp/,${addsuffix /b.ps,${BASES}}}

# # ----------------------------------------------------------------------
# # Used for the gl/harm/d006 basis so that it fits better in the page. Excludes
# # element 0 (the unit constant) leaving 48 instead of 49 elements.
# 
# ELEMS := \
#   001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 \
#   018 019 020 021 022 023 024 025 026 027 028 029 030 031 032 033 034 \
#   035 036 037 038 039 040 041 042 043 044 045 046 047 048
# 
# PSOPTIONS :=  \
#   -ps -paperSize a3 \
#   -figSize 35 -meshSize 0.8 \
#   -obs 1 10 10 10 \
#   -fRange 1.0 -fStep 0.2
# 
# EPSOPTIONS :=  \
#   -eps \
#   -figSize 50 -meshSize 1.0 \
#   -obs 1 10 10 10 \
#   -fRange 1.0 -fStep 0.2
# 
# SELECT := ${addprefix -index ,${ELEMS}}
# # ----------------------------------------------------------------------

# # ----------------------------------------------------------------------
# # Used for the gl/harm/d017 basis
# 
# ELEMS := \
#   000 001 002 003 004 005 006 \
#   025 030 045 050 055 060 \
#   120 130 140 150 \
#   200 250 275 \
#   320 321 322 323
# 
# PSOPTIONS :=  \
#   -ps -paperSize a3 \
#   -figSize 35 -meshSize 0.8 \
#   -obs 1 10 10 10 \
#   -fRange 1.0 -fStep 0.2 \
#   -caption '%I'
# 
# EPSOPTIONS :=  \
#   -eps \
#   -figSize 50 -meshSize 1.0 \
#   -obs 1 10 10 10 \
#   -fRange 1.0 -fStep 0.2
# 
# SELECT := ${addprefix -index ,${ELEMS}}
# # ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Used for a typical spline basis. Plots al elements connected with 
# face 0, edge 0, and vertex 0.

SELECT := \
  -edge 0 -vertex 0 -face 0

PSOPTIONS := \
  -showTriang \
  -ps -paperSize letter \
  -figSize 35 -meshSize 0.600 \
  -radius 1.3 \
  -caption '%I' \
  -caption "step = %D"

EPSOPTIONS := \
  -showTriang \
  -eps \
  -figSize 50 -meshSize 0.750 \
  -radius 1.3
# ----------------------------------------------------------------------

all: ${PSFILES} 

# PSVIEW := kghostview
PSVIEW := evince

${DATADIR}/sp/%/b.ps: ${DATADIR}/sp/%/b.bas ${PROGDIR}/SPPlotBasis 
	${PROGDIR}/SPPlotBasis \
          -basisName ${DATADIR}/sp/$*/b \
          ${SELECT} \
          -outName ${DATADIR}/sp/$*/b \
          ${PSOPTIONS}
	${PROGDIR}/SPPlotBasis \
          -basisName ${DATADIR}/sp/$*/b \
          ${SELECT} \
          -outName ${DATADIR}/sp/$*/b \
          ${EPSOPTIONS}
	${PSVIEW} ${DATADIR}/sp/$*/b.ps


# Rebuild the program if needed:

${PROGDIR}/SPPlotBasis: ${PROGDIR}
	( cd ${PROGDIR} && make SPPlotBasis )
