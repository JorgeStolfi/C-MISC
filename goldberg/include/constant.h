// -*- c++ -*-
// $Id: constant.h,v 1.3 1997/08/15 17:51:10 pmartin Exp $

#ifndef CONSTANT_H
#define CONSTANT_H
#endif

#define ERROR_CHECK


#define TRUE  1
#define FALSE 0 
#define nil 0
#define VERY_FAR          2147483647
#define LARGE             (VERY_FAR/2)
#define LARGE_FLOW        VERY_FAR
#define MAX_CAPACITY      VERY_FAR

#define NOT_IN_TREE -1

/* status of node regarding to queue */ 
#define OUT_OF_QUEUE   0
#define INACTIVE       1
#define ACTIVE         2

// Parser Constants

#ifndef OFFSET
#define OFFSET 0
#endif

#define MAXLINE       100	/* max line length in the input file */

#if defined SHORTEST_PATH
  #define QUEUE
  #define ARC_FIELDS      3	/* no of fields in arc input  */
  #define P_FIELDS        3       /* no of fields in problem line */
  #define PROBLEM_TYPE "sp"       /* name of problem type*/
  #define MIN_NODE_LINES 1
  #define MAX_NODE_LINES 1
#endif

