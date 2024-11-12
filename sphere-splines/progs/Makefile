# Last edited on 2011-12-12 12:00:40 by stolfi

PROGS := \
  SPPlotFunction \
  SPMakeProcFunction \
  SPMakeHarmonicFunction \
  \
  SPComputeApprox \
  \
  SPMakeTriang \
  SPMake3Triang \
  SPMake57Triang \
  SPPlotTriang \
  SPMakePolyBasis \
  SPMakeHarmonicBasis \
  SPFindBezSplineBasis \
  SPOrthoBasis \
  SPPlotBasis \
  SPExtractBasis \
  SPComputeMatrices \
  SPComputeErrorMap \
  SPUniSolve \
  SPMultiSolve \
   \
  SPSplineTest \
  SPIntTest \
  SPDotTest \
  SPHessTest \
  SPGradTest \
  SPCholeskyTest \
  SPCompareMatrices

PROGS_HOLD_FOR_A_SECOND :=

PROGS_HOLD :=

PROGS_NOT_YET := \
  SPTimeOnlySolveGen \
  SPTimeSpaceSolveGen \
  SPTimeSpaceSolveG1C0 \
  SPTimeOnlySolveG3C1 \
  SPPlotTimeSpaceFunction \
  SPMakeTimeOnlyProcFunction \
  SPMakeTimeSpaceProcFunction \
  SPMatrixDecompTest

LIB_HFILES_NOT_YET :=

LIB_CFILES_NOT_YET :=

IGNORE_PROGS := \
  ${PROGS_HOLD_FOR_A_SECOND} \
  ${PROGS_HOLD} \
  ${PROGS_NOT_YET}

IGNORE := \
  ${addsufix .c, ${IGNORE_PROGS}} \
  ${LIB_HFILES_NOT_YET} \
  ${LIB_CFILES_NOT_YET}

include ${STOLFIHOME}/programs/c/GENERIC-PROGS.make
