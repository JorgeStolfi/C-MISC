// Hey EMACS! This is -*- c++ -*- 
// $Id :$ 

#ifndef STATISTICS_H
#define STATISTICS_H 

class statistics
{
public:
  long relabelCnt;
  long pushCnt;
  long mVal;
  long updateCnt;
  long gNodeCnt;
  long greedCnt;
  long numScans;
  long phaseCnt;        
  long contractIters;   
  long zeroPaths;
  long zeroRelabelCnt;
  long linkCnt;
  long cutCnt;
  statistics()  {relabelCnt=0;pushCnt=0;mVal=0;updateCnt=0;gNodeCnt=0;
                 greedCnt=0; numScans=0; phaseCnt=0; contractIters=0;
		 zeroPaths=0; linkCnt=0; cutCnt=0; zeroRelabelCnt=0;}
};

#ifdef MAIN
statistics stats;
#else
extern statistics stats;
#endif

#endif // End of statistics.h
