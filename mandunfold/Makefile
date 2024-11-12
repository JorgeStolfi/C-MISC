# Last edited on 2011-06-06 21:15:14 by stolfi

PROG := mandunfold

JS_LIBS := \
  libjs.a

all: build run

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make

.PHONY:: single-test

# Limit of ratio of radii of successive {M#0} disks on the real axis:
FEIGENBAUM_RATIO := 4.669201609102990

# Western limit of {M#0} disk centers on the the real axis:
W_TIP_RE := -1.40115519

# Northern limit of {M#0} disk centers on the imag axis for {h(c) = c^2 + 1/4},
# namely {sqrt(1/4 + ${W_TIP_RE})}:
SQRT_W_TIP_IM := 1.2849728

# Whole Mandelbrot:
M_CENTER := 00.00000000 00.00000000
M_RADIUS := 02.01000000 02.01000000

# # Western tip:
# M_CENTER := ${W_TIP_RE} 00.00000000
# M_RADIUS := 00.00000046692 00.00000046692

TESTS := \
  h6-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:01000:6: \
  h7-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:01000:7:
  
OTHER_TESTS := \
  h0-Mtot-A:00.00000000!00.00000000:02.01000000!02.01000000:1024:00300:0: \
  h0-Mcsp-A:+0.25000000!00.00000000:00.01000000!00.01000000:1024:00300:0: \
  h0-Mtip-F:${W_TIP_RE}!00.00000000:00.00001000!00.00001000:1024:10000:0: \
  h0-Mtip-H:${W_TIP_RE}!00.00000000:00.00000010!00.00000010:1024:90000:0: \
   \
  h0-M01a-C:-1.75000000!00.00000000:00.05000000!00.03000000:1024:00300:0: \
  h0-M01b-C:-1.75000000!00.00000000:03.00000000!03.00000000:1024:00300:0: \
   \
  h1-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:00300:1: \
  h2-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:00300:2: \
  h3-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:00300:3: \
  h4-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:00300:4: \
  h5-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:10000:5: \
  h6-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:00300:6: \
   \
  h7-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:00300:7: \
  h7-Wtip-E:+3.14592654!00.00000000:00.00010000!00.00010000:0256:50000:7: \
   \
  h8-Wtot-Z:00.00000000!00.00000000:00.00000000!00.00000000:1024:00300:8: \
  
run: all-tests

all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt single-test ; \
        done
        
########################################################################
# For recursive make - user must define ${TEST}

ifneq "/${TEST}" "/"

BANG := !

TARGS := ${subst :, ,${TEST}}

TAG :=      ${word 1,${TARGS}}
W_CENTER := ${subst ${BANG}, ,${word 2, ${TARGS}}}
W_RADIUS := ${subst ${BANG}, ,${word 3, ${TARGS}}}
SIZE :=     ${word 4, ${TARGS}}
NITER :=    ${word 5, ${TARGS}}
FUNC :=     ${word 6, ${TARGS}}

OUT_DIR := out
OUT_NAME := ${OUT_DIR}/${TAG}

OUT_NAME_W := ${OUT_NAME}-W
PGM_FILE_W := ${addsuffix .pgm, ${OUT_NAME_W}}
PNG_FILE_W := ${addsuffix .png, ${OUT_NAME_W}}

OUT_NAME_M := ${OUT_NAME}-M
PGM_FILE_M := ${addsuffix .pgm, ${OUT_NAME_M}}
PNG_FILE_M := ${addsuffix .png, ${OUT_NAME_M}}

single-test: ${PNG_FILE_W}

${PNG_FILE_W}: ${PROG} Makefile
	mkdir -p ${OUT_DIR}
	${PROG} \
            -func ${FUNC} \
            -size ${SIZE} \
            -domain ${W_CENTER} ${W_RADIUS} \
            -range  ${M_CENTER} ${M_RADIUS} \
            -niter ${NITER} \
            -outName ${OUT_NAME}
	cat ${PGM_FILE_W} | pnmtopng  > ${PNG_FILE_W}
	cat ${PGM_FILE_M} | pnmtopng  > ${PNG_FILE_M}
	display ${PNG_FILE_M} & display ${PNG_FILE_W}

endif
# End of ${TEST} section
########################################################################
