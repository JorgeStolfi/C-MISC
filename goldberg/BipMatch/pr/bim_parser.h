// -*- c++ -*-

#ifndef BIM_PARSER_H
#define BIM_PARSER_H

/***************************************************************************/
/*                                                                         */
/*            Parser functions for bipartite matching problems             */
/*                                                                         */
/***************************************************************************/

#include <stdio.h>
#include <string.h>

#define REVERSE_ARCS

#ifdef  NO_ERROR_HANDLING
#define throw(string)   {cout << string;exit(1);}
#define try
#endif

#define ARTIFICIAL_SINK
#define ARTIFICIAL_SOURCE
#define ARC_FIELDS      3	      // no of fields in arc input 
#define NODE_FIELDS     2             // no of fields in node input
#define P_FIELDS        3             // no of fields in problem line 
#define PROBLEM_TYPE "max"            // name of problem type
#define MIN_NODE_LINES 0
#define MAX_NODE_LINES (G.numArcs())
#undef  OFFSET
#define OFFSET 1

static bool_t flowInput;
static long XSize;

inline void swap(long &head,long &tail) 
{
  long tmp;
  tmp =head;
  head = tail;
  tail = tmp;
}

bool_t reverse = false;
inline void readProblemLine(graph &G, char *inputLine)
{
  char    prType[20];     // For reading type of the problem    
  node    *v;
  long    a1,a2,a3,m,n;
  int     fields;

  fields = sscanf ( inputLine, "%*c %20s %ld %ld %ld", prType, &a1, &a2, &a3 ); 
//  n+=2;
  if ( !strcmp ( prType, "max" ) ) {
    flowInput = true;
    if (fields != 3)
      throw("wrong number of parameters in the problem line");
    n = a1;
    m = a2;
  }
  else if (!strncmp ( prType, "bipartite-matching",9) ) {
    flowInput = false;
    if (fields != 4)
      throw("wrong number of parameters in the problem line");
    n = a1+a2+2;
    XSize = a1;
    if (reverse) XSize=a2;
    m = a3+a1+a2;
  }
  else
    throw("It is not a Bipartite Matching problem line");

#ifdef REVERSE_ARCS
  m = m*2;
#endif

  if ( n <= 0  || m <= 0 )
    throw("bad value of a parameter in the problem line");
    
// Allocating Memory For  'nodes', 'arcs'
  G.makeNodes(n);
  G.makeArcs(m);
  if (!flowInput) {
    G.setSource(G.node_i(n-2));
    G.setSink(G.node_i(n-1));
  }
}

inline void readNodeDescription(graph &G, char *input_line)
{
  if (!flowInput)
    throw("Node line not allowed in bipartite graph formulation");

  long i;
  char c;

  // Reading Node Description
  if (sscanf ( input_line,"%*c %ld %c", &i, &c ) < NODE_FIELDS)
    throw("wrong number of parameters in a node description line");
  i -= OFFSET;
  
  if ( i < 1 || i >= G.numNodes() )
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
  int fields;

  fields = sscanf ( input_line,"%*c %ld %ld %ld", &tail, &head, &ub);
  if (fields != ((flowInput) ? 3 : 2))
    throw("wrong number of parameters in an arc desrciption line");
  tail -= OFFSET;
  head -= OFFSET;

  if (reverse) swap(head,tail);
  if (!flowInput) head += XSize;
  //cout << head << " " << tail << "\n";

  if ( tail < 0  ||  tail >= G.numNodes() || 
       head < 0  ||  head >= G.numNodes() )
    throw("invalid node reference in an arc description line");

//  if ( ub < 0) ub = BIGGEST_FLOW;
  if (flowInput && ub != 1) cout << "Warning Upper Bound changed to 1\n";

  e->setHead(G.node_i(head));
  (e+1)->setHead(G.node_i(tail));

  return tail;
}

inline bool_t connectedToSource(node *v)
{
  if (flowInput || (v->index() >= XSize)) return false;
  else return true;
}

inline bool_t connectedToSink(node *v)
{
  if (flowInput || (v->index() < XSize) || (v->index() > G.numNodes()-3)) return false;
  else return true;
}

inline void createSourceArc(arc *e,node *v)
{
  //cout << G.getSource()->index() << " " << v->index() << "\n";
  e->setHead(v);
  (e+1)->setHead(G.getSource());
}

inline void createSinkArc(arc *e,node *v)
{
  //  cout << G.getSink()->index() << " " << v->index() << "\n";
  e->setHead(G.getSink());
  (e+1)->setHead(v);
}

#endif // End of bim_parser.h

