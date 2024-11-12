# Last edited on 2015-05-04 15:23:23 by stolfilocal

PROG := fit_btc_bubble_model
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

#all: run-prog1
all: run-prog2

DATA_DIR := data
OUT_DIR := out
# HRAD := 07
HRAD := 00

# RUN_DATE := 2015-04-28
# NB := 10

RUN_DATE := 2015-04-28
NB := 11

# RUN_DATE := 2015-04-28
# NB := 12

# RUN_DATE := 2015-04-28
# NB = 19

# RUN_DATE := 2015-04-28
# NB := 20

# RUN_DATE := 2015-04-28
# NB := 21

# RUN_DATE := 2015-04-23
# NB := 19

# EVAL_RANGE_OPT := 
# EVAL_RANGE_OPT := -evalRange 2010-07-20 2010-12-31
# EVAL_RANGE_OPT := -evalRange 2010-07-20 2011-06-30
# EVAL_RANGE_OPT := -evalRange 2010-11-01 2011-02-18
# EVAL_RANGE_OPT := -evalRange 2010-11-01 2011-06-30
# EVAL_RANGE_OPT := -evalRange 2010-11-01 2012-04-30
# EVAL_RANGE_OPT := -evalRange 2011-01-01 2011-06-30
# EVAL_RANGE_OPT := -evalRange 2011-04-01 2012-03-31
# EVAL_RANGE_OPT := -evalRange 2011-11-01 2012-10-30
# EVAL_RANGE_OPT := -evalRange 2012-07-01 2013-07-31
# EVAL_RANGE_OPT := -evalRange 2013-01-01 2013-11-30
# EVAL_RANGE_OPT := -evalRange 2013-01-01 2015-02-18
# EVAL_RANGE_OPT := -evalRange 2013-02-01 2013-11-30
# EVAL_RANGE_OPT := -evalRange 2013-04-01 2013-12-31
EVAL_RANGE_OPT := -evalRange 2013-09-01 2015-04-13
# EVAL_RANGE_OPT := -evalRange 2013-11-01 2014-07-31
# EVAL_RANGE_OPT := -evalRange 2013-11-15 2015-02-15
# EVAL_RANGE_OPT := -evalRange 2014-02-01 2014-07-31
# EVAL_RANGE_OPT := -evalRange 2014-04-01 2014-06-30
# EVAL_RANGE_OPT := -evalRange 2014-09-01 2014-12-31
# EVAL_RANGE_OPT := -evalRange 2014-09-01 2015-03-15
# EVAL_RANGE_OPT := -evalRange 2014-12-01 2015-04-13

IN_PRICE_FILE := ${DATA_DIR}/2010-07-17--2015-04-13-PREF-USD-01d-sm${HRAD}.txt

IN_PARMS_PREFIX := ${DATA_DIR}/${RUN_DATE}-bubble-${NB}
IN_PARMS_FILE := ${IN_PARMS_PREFIX}.parms
IN_LO_PARMS_FILE := ${IN_PARMS_PREFIX}-lo.parms
IN_HI_PARMS_FILE := ${IN_PARMS_PREFIX}-hi.parms

OUT_PREFIX_1 := ${OUT_DIR}/${RUN_DATE}-bubble-${NB}-sm${HRAD}-run1
OUT_PREFIX_2 := ${OUT_DIR}/${RUN_DATE}-bubble-${NB}-sm${HRAD}-run2
OUT_SUFFIX := fit


OUT_PRICE_FILE_1 := ${OUT_PREFIX_1}-${OUT_SUFFIX}.txt

run-prog1: ${PROGDIR}/${PROG} ${OUT_DIR} ${IN_PRICE_FILE} ${IN_PARMS_FILE}
	${PROGDIR}/${PROG} \
            -inParms ${IN_PARMS_FILE} \
            -maxLSQIters 2 \
            -hrad ${HRAD} \
            -truncate 0.050 \
            -outPrefix ${OUT_PREFIX_1} \
          < ${IN_PRICE_FILE} \
          > ${OUT_PRICE_FILE_1}
	plot_single_bubbles.sh ${RUN_DATE} ${HRAD} ${NB} ${OUT_PRICE_FILE_1}
	plot_bubbles.sh ${RUN_DATE} ${HRAD} ${NB} ${OUT_PRICE_FILE_1}

OUT_PRICE_FILE_2 := ${OUT_PREFIX_2}-${OUT_SUFFIX}.txt

run-prog2: ${PROGDIR}/${PROG} ${OUT_DIR} ${IN_PRICE_FILE} ${IN_PARMS_FILE}
	${PROGDIR}/${PROG} \
            -inParms ${IN_PARMS_FILE} \
            -inRangeParms ${IN_LO_PARMS_FILE} ${IN_HI_PARMS_FILE} \
            ${EVAL_RANGE_OPT} \
            -maxLSQIters 4 \
            -maxNLIters 5 \
            -hrad ${HRAD} \
            -truncate 0.050 \
            -outPrefix ${OUT_PREFIX_2} \
          < ${IN_PRICE_FILE} \
          > ${OUT_PRICE_FILE_2}
	plot_single_bubbles.sh ${RUN_DATE} ${HRAD} ${NB} ${OUT_PRICE_FILE_2}
	plot_bubbles.sh ${RUN_DATE} ${HRAD} ${NB} ${OUT_PRICE_FILE_2}
          
          

clean::

