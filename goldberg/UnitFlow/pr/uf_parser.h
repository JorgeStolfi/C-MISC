// -*- c++ -*-

#ifndef UF_PARSER_H
#define UF_PARSER_H

/***************************************************************************/
/*                                                                         */
/*            Parser functions for unit-capacity flow problems             */
/*                                                                         */
/***************************************************************************/

#include <stdio.h>

#ifdef  NO_ERROR_HANDLING
#ifndef throw
#define throw(string)  {cout << string<<"\n";exit(1);}  
#define try
#endif
#endif

#define ARC_FIELDS      3	      // no of fields in arc input 
#define NODE_FIELDS     2             // no of fields in node input
#define P_FIELDS        3             // no of fields in problem line 
#define PROBLEM_TYPE "max"            // name of problem type
#define PROBLEM_TYPE2 "unit-flow"      // name of problem type
#define MIN_NODE_LINES 1
#define MAX_NODE_LINES (G.numArcs())

inline void readProblemLine(graph &G, char *inputLine)
{
  char    prType[10];     // For reading type of the problem    
  node    *v;
  long    m,n;

  if (sscanf ( inputLine, "%*c %s %ld %ld", prType, &n, &m ) != P_FIELDS) 
    throw("wrong number of parameters in the problem line");
  if ( strcmp ( prType, PROBLEM_TYPE ) && strcmp ( prType, PROBLEM_TYPE2) ) 
    throw("it is not a Unit Flow problem line");

#ifdef REVERSE_ARCS
  m = m*2;
#endif

  if ( n <= 0  || m <= 0 )
    throw("bad value of a parameter in the problem line");
    
// Allocating Memory For  'nodes', 'arcs'
  G.makeNodes(n);
  G.makeArcs(m);
}

inline void readNodeDescription(graph &G, char *input_line)
{
  long i;
  char c;

  // Reading Node Description
  if (sscanf ( input_line,"%*c %ld %c", &i, &c ) < NODE_FIELDS)
    throw("wrong number of parameters in a node description line");
  i -= param.offset;
  
  if ( i < 0 || i >= G.numNodes() )
    throw("invalid node number");
  
  switch(c) 
    {
    case 't':
      if (G.getSink() != nil) throw("Multiple sinks not permitted");
      G.setSink(G.node_i(i));
      break;

    case 's':
      if (G.getSource() != nil) throw("Multiple sources not permitted");
      G.setSource(G.node_i(i));
      break;

    default:
      throw("Improper source-sink description in node description line");
    }
}


inline long readArcDescription(graph &G, char *input_line, arc *e)
{
  long head, tail, ub;
  
  if ( sscanf ( input_line,"%*c %ld %ld %ld", &tail, &head, &ub)
      != ARC_FIELDS )
        throw("wrong number of parameters in an arc desrciption line");
  tail -= param.offset;
  head -= param.offset;
		
  if ( tail < 0  ||  tail >= G.numNodes() || 
       head < 0  ||  head >= G.numNodes() )
    throw("invalid node reference in an arc description line");

//  if ( ub < 0) ub = BIGGEST_FLOW;
//  if ( ub != 1) cout << "Warning Upper Bound changed to 1\n";

  e->setHead(G.node_i(head));
  e->setReverse(e+1);
  e->setUsed(false);
  e->setForward(true);
  (e+1)->setHead(G.node_i(tail));
  (e+1)->setReverse(e);
  (e+1)->setUsed(true);
  (e+1)->setForward(false);
  return tail;
}

#endif UF_PARSER_H
