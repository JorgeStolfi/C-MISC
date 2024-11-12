# Last edited on 2023-02-04 07:03:11 by stolfi

PROGS := \
  SPFindBezSplineBasis \
  SPOrthoBasis \
  SPPlotBasis
  
PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

SMPORDER := 40
SO := ${SMPORDER}

# TRI:HTAG:DEG:CONT:STYLE:ORTHO:SHOW

RUNS := \
  icosa/r57/zy/01:h:2:-1:NEW:YES:f0,e0,v0 \
  icosa/r57/zy/01:nh:1:-1:NEW:YES:f0,e0,v0 \
   \
  icosa/r57/zy/01:h:3:-1:NEW:YES:f0,e0,v0 \
  icosa/r57/zy/01:nh:2:-1:NEW:YES:f0,e0,v0

RUNS_DISCONTINUOUS := \
  icosa/r57/zy/01:h:2:-1:NEW:NO:f0,e0,v0 \
  icosa/r57/zy/01:nh:1:-1:NEW:NO:f0,e0,v0 \
   \
  icosa/r57/zy/01:h:3:-1:NEW:NO:f0,e0,v0 \
  icosa/r57/zy/01:nh:2:-1:NEW:NO:f0,e0,v0

RUNS_MAYBE := \
  icosa/raw:h:8:1:OLD:NO:f0,e0,v0 \
  icosa/raw:h:8:1:OLD:YES:f0,e0,v0

RUNS_57 := \
  icosa/r57/zy/01:h:7:1:NEW:NO:f0,e0,v0 \
  icosa/r57/zy/01:h:7:1:NEW:YES:f0,e0,v0

DONE := \
  icosa/raw:nh:7:1:NEW:NO:f0,e0,v0 \
  icosa/raw:nh:7:1:NEW:YES:f0,e0,v0 \
   \
  icosa/r57/zz/01:nh:6:1:NEW:NO:f0,e0,v0 \
  icosa/r57/zz/01:nh:6:1:NEW:YES:f0,e0,v0 \
   \
  icosa/r57/zy/01:nh:6:1:NEW:NO:f0,e0,v0 \
  icosa/r57/zy/01:nh:6:1:NEW:YES:f0,e0,v0 \
   \
  icosa/raw:nhx:6:1:NEW:YES:f0,e0,v0 \
   \
  icosa/raw:h:7:1:OLD:NO:f0,e0,v0 \
  icosa/raw:h:7:1:OLD:YES:f0,e0,v0 \
  icosa/raw:h:5:1:NEW:YES:f0,e0,v0 \
  icosa/raw:h:5:1:NEW:NO:f0,e0,v0 \
  icosa/raw:nh:6:1:NEW:NO:f0,e0,v0 \
  icosa/raw:nh:6:1:NEW:YES:f0,e0,v0

NOTESTRUNS := \
  icosa/raw:h:7:1:OLD:NO:f0,e0,v0 \
  icosa/raw:h:7:1:OLD:YES:f0,e0,v0 \
   \
  icosa/raw:h:7:1:OLD:NO:f0,e0,v0 \
  icosa/raw:nh:6:1:OLD:NO:f0,e0,v0 \
  \
  icosa/raw:h:5:1:NEW:YES:f0,e0,v0 \
  icosa/raw:h:6:1:NEW:YES:f0,e0,v0 \
  icosa/raw:h:7:1:NEW:YES:f0,e0,v0 \
  icosa/raw:nh:6:1:NEW:YES:f0,e0,v0

RUNSRAW := \
  icosa/raw:h:0:-1:NEW:YES:f0,e0,v0 \
  icosa/raw:h:1:-1:NEW:YES:f0,e0,v0 \
   \
  icosa/raw:h:1:0:NEW:YES:f0,e0,v0 \
  icosa/raw:h:2:0:NEW:YES:f0,e0,v0 \
  icosa/raw:h:3:0:NEW:YES:f0,e0,v0 \
  icosa/raw:h:4:0:NEW:YES:f0,e0,v0 \
  icosa/raw:h:5:0:NEW:YES:f0,e0,v0 \
   \
  icosa/raw:h:5:1:NEW:YES:f0,e0,v0 \
  icosa/raw:h:6:1:NEW:YES:f0,e0,v0 \
  icosa/raw:h:7:1:NEW:YES:f0,e0,v0 \
   \
  icosa/raw:nh:1:-1:NEW:YES:f0,e0,v0 \
   \
  icosa/raw:nh:2:0:NEW:YES:f0,e0,v0 \
  icosa/raw:nh:3:0:NEW:YES:f0,e0,v0 \
  icosa/raw:nh:4:0:NEW:YES:f0,e0,v0 \
  icosa/raw:nh:5:0:NEW:YES:f0,e0,v0 \
   \
  icosa/raw:nh:6:1:NEW:YES:f0,e0,v0 \
  icosa/raw:nh:7:1:NEW:YES:f0,e0,v0
  
RUNSREF := \
   \
  icosa/r3/01:h:4:0:NEW:YES:f0,e0,e30,v0,v20 \
  icosa/r3/01:nh:3:0:NEW:YES:f0,e0,e30,v0,v20 \
  \
  icosa/r3/01:h:7:1:NEW:YES:f0,e0,e30,v0,v20 \
  icosa/r3/01:nh:6:1:NEW:YES:f0,e0,e30,v0,v20 \
  \
  icosa/r3/02:h:7:1:NEW:YES:f0,e0,e30,v0,v20 \
  icosa/r3/02:nh:6:1:NEW:YES:f0,e0,e30,v0,v20 \
  \
  icosa/r3/03:h:7:1:NEW:YES:0,e0,e30,v0,v20 \
  icosa/r3/03:nh:6:1:NEW:YES:f0,e0,e30,v0,v20 \

RUNSBUG := \
  icosa/raw:nh:1:-1:NEW:YES:f0,e0,v0 \
  icosa/r3/03:nh:6:1:NEW:YES:f0,e0,e30,v0,v20

# RANGEOPTS := -fRange 2.5 -fStep 0.5
# RANGEOPTS := -fRange 0.5 -fStep 0.10 
# RANGEOPTS := -fRange 1.5 -fStep 0.3  -eps
  
.PHONY:: \
  clean-all run-all \
  clean-single run-single \
  make-basis plot-basis show-elems show-elem

all: clean-all run-all
# all: run-all
# all: clean-all

clean-all: 
	for r in ${RUNS}; do \
          ${MAKE} RUNARGS="$$r" clean-single ; \
        done

run-all: 
	for r in ${RUNS}; do \
          ${MAKE} RUNARGS="$$r" run-single ; \
        done

# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS := ${subst :, ,${RUNARGS}}
TRI   := ${word 1,${FIELDS}}
HTAG  := ${word 2,${FIELDS}}
DEG   := ${word 3,${FIELDS}}
CONT  := ${word 4,${FIELDS}}
STYLE := ${word 5,${FIELDS}}
ORTHO := ${word 6,${FIELDS}}
SHOW  := ${word 7,${FIELDS}}

CONTTAG := ${subst -1,n,${CONT}}

# Should we use the old-style ANS basis?
ifeq "/${STYLE}" "/NEW"
  STYLETAG := 
  STYLEOPT := -newStyle
else
  STYLETAG := -OLD
  STYLEOPT := -oldStyle
endif

# Plotting options:
COMMA := ,
SHOW := ${subst ${COMMA}, @,@${SHOW}}
SHOW := ${subst @f, -face ,${SHOW}}
SHOW := ${subst @e, -edge ,${SHOW}}
SHOW := ${subst @v, -vertex ,${SHOW}}
SHOW := ${subst @i, -index ,${SHOW}}

TRIDIR := ${DATADIR}/sp/pw/${TRI}
TRINAME := ${TRIDIR}/t

# Should we orthogonalize the basis?
ifeq "/${ORTHO}" "/NO"
  ORTHOTAG := -NO
else
  ORTHOTAG := 
endif

OUTDIR := ${TRIDIR}/${HTAG}-d${DEG}c${CONTTAG}${STYLETAG}${ORTHOTAG}
OUTNAME := ${OUTDIR}/b

clean-single: 
	-/bin/rm -f ${OUTNAME}.bas ${OUTNAME}-x.bas
	-/bin/rm -f ${OUTNAME}-el0*.sfn ${OUTNAME}-el0*.eps ${OUTNAME}.ps

run-single: clean-single make-basis plot-basis show-elems
# run-single: clean-single make-basis plot-basis

make-basis: ${OUTNAME}.bas

${OUTDIR}:
	mkdir -p ${OUTDIR}

ifeq "/${ORTHO}" "/NO"
  ${OUTNAME}.bas: ${PROGDIR}/SPFindBezSplineBasis Makefile \
            ${OUTDIR} ${TRINAME}.tri
	${PROGDIR}/SPFindBezSplineBasis \
            -triName ${TRINAME} \
            ${subst h,-homogeneous,${subst nh,,${HTAG}}} \
            -degree ${DEG} -continuity ${CONT} ${STYLEOPT} \
            -outName ${OUTNAME} \
          2>&1 | tee ${OUTNAME}.log
else
  INDIR := ${TRIDIR}/${HTAG}-d${DEG}c${CONTTAG}${STYLETAG}-NO
  INNAME := ${INDIR}/b

  # ORTHOOPTS := -sortGDist -gram 1 -eigen 1 
  ORTHOOPTS := -sortGDist -gram 1 -eigen 1 -precise
  # ORTHOOPTS := -eigen 2
  # ORTHOOPTS := -eigen 3 -gram 1
  # ORTHOOPTS := -gram 2

  ${OUTNAME}.bas: ${PROGDIR}/SPOrthoBasis Makefile \
            ${INNAME}.bas ${OUTDIR} ${TRINAME}.tri
	${PROGDIR}/SPOrthoBasis \
            -basisName ${INNAME} \
            ${ORTHOOPTS} \
            -smpOrder ${SMPORDER} \
            -checkOrtho 1.0e-10 \
            -outName ${OUTNAME} \
          2>&1 | tee ${OUTNAME}.log
endif

plot-basis: ${PROGDIR}/SPPlotBasis \
            ${OUTNAME}.bas Makefile
	${PROGDIR}/SPPlotBasis \
          -basisName ${OUTNAME} \
          ${SHOW} \
          ${RANGEOPTS} \
          -autoObs 0.90 -showTriang -showSupp \
          -meshSize 1.5 -figSize 40 -paperSize a3 \
          -outName ${OUTNAME}

# PSVIEW := kghostview
PSVIEW := evince

show-elems: 
	for el in `ls ${OUTNAME}.ps ${OUTNAME}-el0*.eps` ; do \
          ${PSVIEW} $$el ; \
        done

endif
# End of ${RUNARGS} part.
# ======================================================================
