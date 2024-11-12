# Last edited on 2013-10-26 18:51:47 by stolfilocal

PROG := compute_interpolators
PROGDIR := ..

OUT_PREFIX := out/splines
GPL_FILE := ${OUT_PREFIX}.gpl

MAKEPLOTS := makeplots-B
# MAKEPLOTS := makeplots-I
# MAKEPLOTS := makeplots-O

run: ${GPL_FILE} ${MAKEPLOTS}

${GPL_FILE}: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} > ${GPL_FILE}
	cat ${GPL_FILE}

makeplots-B: ${GPL_FILE}
	plot-spline.sh ${GPL_FILE} "D0" Bn -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D0" B0 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" B0 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D0" B1 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" B1 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D2" B1 -2.5 +2.5
	plot-spline.sh ${GPL_FILE} "D0" B2 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" B2 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D2" B2 -2.5 +2.5
	plot-spline.sh ${GPL_FILE} "D3" B2 -3.5 +3.5

makeplots-I: ${GPL_FILE}
	plot-spline.sh ${GPL_FILE} "D0" In -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D0" I0 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" I0 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D0" I1 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" I1 -2.0 +2.0
	plot-spline.sh ${GPL_FILE} "D2" I1 -4.5 +4.5
	plot-spline.sh ${GPL_FILE} "D0" I2 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" I2 -2.0 +2.0
	plot-spline.sh ${GPL_FILE} "D2" I2 -5.5 +5.5
	plot-spline.sh ${GPL_FILE} "D3" I2 -9.9 +9.9

makeplots-O: ${GPL_FILE}
	plot-spline.sh ${GPL_FILE} "D0" O0 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" O0 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D0" O1 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" O1 -2.0 +2.0
	plot-spline.sh ${GPL_FILE} "D2" O1 -4.5 +4.5
	plot-spline.sh ${GPL_FILE} "D0" O2 -1.5 +1.5
	plot-spline.sh ${GPL_FILE} "D1" O2 -2.0 +2.0
	plot-spline.sh ${GPL_FILE} "D2" O2 -5.5 +5.5
	plot-spline.sh ${GPL_FILE} "D3" O2 -9.9 +9.9
        
include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
