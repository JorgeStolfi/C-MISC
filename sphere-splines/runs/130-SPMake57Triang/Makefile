# Last edited on 2023-02-04 07:03:52 by stolfi

PROGS := \
  SPMake57Triang \
  SPPlotTriang
  
PROGDIR := ../../progs

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DATADIR := ../../data

# TRI:TAGIN:TAGOT:TYPE:DEPTH:NOPTS:WAREA:WLENGTH:WANGLE:WXRATIO:AVGMAX

RUNS := \
  icosa:raw:r57:zy:01:00000:1.0000:1.0000:1.0000:0.0000:-average

NONRUNS := \
  icosa:raw:r57:zz:03:00000:1.0000:1.0000:1.0000:0.0000:-average \
  icosa:raw:r57:zz:02:00000:1.0000:1.0000:1.0000:0.0000:-average \
  icosa:raw:r57:zz:01:00000:1.0000:1.0000:1.0000:0.0000:-average \
   \
  icosa:raw:r57:mb:02:10000:0.0000:1.0000:0.0000:0.0000:-maximum \
   \
  icosa:raw:r57:ah:01:10000:1.0000:1.0000:1.0000:0.0000:-average \
  icosa:raw:r57:ah:02:10000:1.0000:1.0000:1.0000:0.0000:-average \
  icosa:raw:r57:ah:03:10000:1.0000:1.0000:1.0000:0.0000:-average \
  \
  icosa:raw:r57:aa:01:10000:1.0000:0.0000:0.0000:0.0000:-average \
  icosa:raw:r57:ma:01:10000:1.0000:0.0000:0.0000:0.0000:-maximum \
   \
  icosa:raw:r57:ab:01:10000:0.0000:1.0000:0.0000:0.0000:-average \
  icosa:raw:r57:mb:01:10000:0.0000:1.0000:0.0000:0.0000:-maximum \
   \
  icosa:raw:r57:ac:01:10000:0.0000:0.0000:1.0000:0.0000:-average \
  icosa:raw:r57:mc:01:10000:0.0000:0.0000:1.0000:0.0000:-maximum \
   \
  icosa:raw:r57:ad:01:10000:0.0000:0.0000:0.0000:1.0000:-average \
  icosa:raw:r57:md:01:10000:0.0000:0.0000:0.0000:1.0000:-maximum \
   \
  icosa:raw:r57:ae:01:10000:1.0000:1.0000:0.0000:0.0000:-average \
  icosa:raw:r57:me:01:10000:1.0000:1.0000:0.0000:0.0000:-maximum \
   \
  icosa:raw:r57:af:01:10000:1.0000:0.0000:1.0000:0.0000:-average \
  icosa:raw:r57:mf:01:10000:1.0000:0.0000:1.0000:0.0000:-maximum \
   \
  icosa:raw:r57:ag:01:10000:1.0000:0.0000:0.0000:1.0000:-average \
  icosa:raw:r57:mg:01:10000:1.0000:0.0000:0.0000:1.0000:-maximum \
   \
  icosa:raw:r57:ah:01:10000:1.0000:1.0000:1.0000:0.0000:-average \
  icosa:raw:r57:mh:01:10000:1.0000:1.0000:1.0000:0.0000:-maximum \
   \
  icosa:raw:r57:ai:01:10000:1.0000:1.0000:0.0000:1.0000:-average \
  icosa:raw:r57:mi:01:10000:1.0000:1.0000:0.0000:1.0000:-maximum \
   \
  icosa:raw:r57:aj:01:10000:1.0000:0.0000:1.0000:1.0000:-average \
  icosa:raw:r57:mj:01:10000:1.0000:0.0000:1.0000:1.0000:-maximum \
   \
  icosa:raw:r57:ak:01:10000:1.0000:1.0000:1.0000:1.0000:-average \
  icosa:raw:r57:mk:01:10000:1.0000:1.0000:1.0000:1.0000:-maximum \

.PHONY:: run clean-single run-single \
  showplots-single showplot-opt showplot-sec

all: run

clean:
	for tri in ${RUNS} ; do \
          ${MAKE}  RUNARGS="$$tri" clean-single ; \
        done

run:
	for tri in ${RUNS} ; do \
          ${MAKE}  RUNARGS="$$tri" run-single ; \
        done

# ======================================================================
# For recursive "make" only -- caller must define ${RUNARGS}

RUNARGS := RUNARGS.IS.UNDEFINED
ifneq "/${RUNARGS}" "/RUNARGS.IS.UNDEFINED"

FIELDS   := ${subst :, ,${RUNARGS}}
TRI      := ${word  1,${FIELDS}}
TAGIN    := ${word  2,${FIELDS}}
TAGOT    := ${word  3,${FIELDS}}
TYPE     := ${word  4,${FIELDS}}
DEPTH    := ${word  5,${FIELDS}}
NOPTS    := ${word  6,${FIELDS}}
WAREA    := ${word  7,${FIELDS}}
WLENGTH  := ${word  8,${FIELDS}}
WANGLE   := ${word  9,${FIELDS}}
WXRATIO  := ${word 10,${FIELDS}}
AVGMAX   := ${word 11,${FIELDS}}

INDIR := ${DATADIR}/sp/pw/${TRI}/${TAGIN}/
INNAME := ${INDIR}/t

OTDIR := ${DATADIR}/sp/pw/${TRI}/${TAGOT}/${TYPE}/${DEPTH}
OTNAME := ${OTDIR}/t
PSNAME := ${OTNAME}-py

clean-single:
	-/bin/rm -f ${OTNAME}.tri ${OTNAME}.log  ${PSNAME}.ps

run-single: clean-single ${OTNAME}.tri ${PSNAME}.ps

${OTDIR}:
	mkdir -p ${OTDIR}

ifeq "${TYPE}:${DEPTH}" "xa:01" 
# Area, depth = 1
SHAPEOPTS := \
  -shape 00 0.3659 0.4615 0.1726   0.6640 0.2741 0.0619
endif

ifeq "${TYPE}:${DEPTH}" "xb:02" 
# Area, depth = 2, maximum, best 1.17305380e-01
SHAPEOPTS := \
  -shape 00  0.3361 0.4807 0.1832  0.6347 0.2683 0.0970 \
  -shape 01  0.3782 0.4582 0.1636  0.6797 0.2191 0.1012
endif

ifeq "${TYPE}:${DEPTH}" "xc:02" 
# Area, depth = 2, average
SHAPEOPTS := \
  -shape 00 0.3558 0.4609 0.1832   0.6312 0.2808 0.0879 \
  -shape 01 0.3670 0.4606 0.1724   0.6681 0.2471 0.0847
endif

ifeq "${TYPE}:${DEPTH}" "xd:02" 
# Length, depth = 2
SHAPEOPTS := \
  -shape 00 0.3558 0.4609 0.1832   0.6312 0.2808 0.0879 \
  -shape 01 0.3670 0.4606 0.1724   0.6681 0.2471 0.0847
endif

ifeq "${TYPE}:${DEPTH}" "xe:02" 
# Area + length + bit of xratio, depth = 2, maximum, best 1.54526829e+00
SHAPEOPTS := \
  -shape 00  0.3536 0.4628 0.1837  0.6779 0.3041 0.0180 \
  -shape 01  0.3701 0.4663 0.1636  0.6777 0.2577 0.0646
endif

ifeq "${TYPE}:${DEPTH}" "zz:03" 
# Fixed, depth = 03
SHAPEOPTS := \
  -shape 00  0.3443 0.4605 0.1952  0.6592 0.2952 0.0456 \
  -shape 01  0.3501 0.4669 0.1830  0.6766 0.2768 0.0466 \
  -shape 02  0.3501 0.4669 0.1830  0.6766 0.2768 0.0466
endif

ifeq "${TYPE}:${DEPTH}" "zy:01" 
# Fixed, depth = 01
SHAPEOPTS := \
  -swapEdges \
  -shape 00  0.3350 0.4650 0.2000  0.6650 0.2900 0.0450
endif

ifeq "${TYPE}:${DEPTH}" "ah:02" 
# ah, depth = 02
SHAPEOPTS := \
  -shape 00  0.3443 0.4605 0.1952  0.6592 0.2952 0.0456 \
  -shape 01  0.3501 0.4669 0.1830  0.6766 0.2768 0.0466
endif

ifeq "${TYPE}:${DEPTH}" "mb:02" 
# mb, depth = 02
SHAPEOPTS := \
  -swapEdges \
  -shape 00  0.3478 0.4448 0.2074  0.6489 0.3141 0.0370 \
  -shape 01  0.3474 0.4698 0.1828  0.6356 0.3044 0.0600
endif

ifeq "${TYPE}:${DEPTH}" "mc:02" 
# # Area + length + bit of xratio, depth = 2, maximum, best 1.10799430e+00
SHAPEOPTS := \
  -swapEdges \
  -shape 00  0.3568 0.4362 0.2070  0.5988 0.2679 0.1333 \
  -shape 01  0.3826 0.4456 0.1718  0.6627 0.2101 0.1272
endif

ifeq "${TYPE}:${DEPTH}" "md:02" 
# Area + length + bit of xratio, depth = 2, maximum, best 4.58719947e-01
SHAPEOPTS := \
  -swapEdges \
  -shape 00  0.3625 0.4574 0.1801  0.7042 0.2910 0.0048 \
  -shape 01  0.3804 0.4618 0.1578  0.6919 0.2409 0.0672
endif

${OTNAME}.tri: ${OTDIR} ${INNAME}.tri \
            ${PROGDIR}/SPMake57Triang Makefile
	${PROGDIR}/SPMake57Triang \
            -triName ${INNAME} \
            -depth ${DEPTH} \
            ${SHAPEOPTS} \
            -optimize ${NOPTS} -minStep 0.00001 -maxStep 0.0100 \
              -area ${WAREA} -length ${WLENGTH} \
              -angle ${WANGLE} -xratio ${WXRATIO} \
              ${AVGMAX} \
            -outName ${OTNAME} \
          2>&1 | tee ${OTNAME}.log

#  -verbose -plot \
# \
#

# PSVIEW := kghostview
PSVIEW := evince

# PAPEROPTIONS := -figSize 230 -paperSize a3
PAPEROPTIONS := -figSize 150 -paperSize letter

${PSNAME}.ps: ${PROGDIR}/SPPlotTriang Makefile \
            ${OTNAME}.tri
	${PROGDIR}/SPPlotTriang \
          -triName ${OTNAME} \
          -outName ${PSNAME} \
          -obs 0 4 2 1 -radius 1.05 \
          -showBack ${PAPEROPTIONS} \
          -caption "${OTNAME}.tri" \
          -caption "${RUNARGS}" \
          -caption "${subst -shape," -caption ",${SHAPEOPTS}}" 
	${PSVIEW} ${PSNAME}.ps

showplots-single: showplot-opt showplot-sec

showplot-opt: ${OUTNAME}-opt.plt m57tri-plot-opt
	( m57tri-plot-opt ${OUTNAME} & )

showplot-sec: ${OUTNAME}-sec.plt m57tri-plot-sec
	( m57tri-plot-sec ${OUTNAME} & )

endif
# End of ${RUNARGS} part.
# ======================================================================
