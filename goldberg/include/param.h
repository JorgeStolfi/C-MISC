// Hey emacs! This is -*- c++ -*-
// $Id: param.h,v 1.3 1997/08/06 19:47:36 mslevine Exp $

#include <constant.h>

#ifndef PARAM_H
#define PARAM_H

class parameters
{
public:
  bool_t    gap;
  double  globalUpdateFreq;
  bool_t    greedy;
  long    offset;
#ifdef P5_COUNTER
  int     p5;
#endif
  parameters() {gap=true; greedy=false; globalUpdateFreq = 1.0; offset=1; 
#ifdef P5_COUNTER
		p5=0;
#endif
	      }
};

#ifdef MAIN
parameters param;
#else
extern parameters param;
#endif

#endif // End of param.h

