// -*- c++ -*-
// $Id: graph.h,v 1.12 1997/07/01 15:25:43 mslevine Exp $
#ifndef GRAPH_H
#define GRAPH_H

/***********************************************************************/
/*                                                                     */
/*         CLASS graph                                                 */
/*                                                                     */
/***********************************************************************/


class graph
{
  node *node_array;      // Start of node array    
  node *src;             // Source node
  node *snk;             // Sink node
  long  n;               // Number of nodes in array $nodes$  
  long n_alloc;          // Number of nodes allocated
  long n_orig;           // original number of nodes
  arc *arc_array;        // Start of arc array
  long  m;               // Number of arcs in array $arcs$
  long m_alloc;          // Number of arcs allocated
  long m_orig;           // original number of arcs
#ifdef FLOW
  CapType maxArcCap;     // maximum arc capcity
#endif

public:

  inline void makeNodes(long num, long alloc);
  inline void makeArcs(long num, long alloc);
  void makeNodes(long num)         { makeNodes(num, num); }
  void makeArcs(long num)          { makeArcs(num, num); }
  long numArcs()                   {return m;}
  long numNodes()                  {return n;}
  inline long name(node *v);
  inline long index(node *v);
  inline long index(arc *e);

  inline node *firstNode();
  inline node *lastNode();
  inline node *succNode(node *v);
  inline node *precNode(node *v);
  inline node *node_i(long i);

  inline node *newNode();

  void freeLastNewNode()           { n--; }
  void freeNewNodes()              { n = n_orig; }

  inline arc *newArc();
  void freeNewArcs()               { m = m_orig; }

  inline arc *firstArc();
  inline arc *lastArc();
  inline arc *succArc(arc *e);
  inline arc *precArc(arc *e);
  inline arc *arc_i(long i);

  inline void setSource(long i);
  inline void setSource(node *v);
  inline node *getSource();

  inline void setSink(long i);
  inline void setSink(node *v);
  inline node *getSink();

#ifdef FLOW
  CapType getMaxArcCapacity()      { return maxArcCap; }
  void setMaxArcCapacity(CapType val) { maxArcCap = val; }

  graph()                          {src=nil; snk=nil; maxArcCap=0;}
#else
  graph()                          {src=nil; snk=nil;}
#endif
  graph(long n_size, long e_size) { makeNodes(n_size);
				    makeArcs(e_size);
				    graph();}

//  void readProblemLine(char *input_line);
//  void readNodeDescription(char *input_line);
//  long readArcDescription(char *input_line, arc *e);
};
#ifdef MAIN
graph G;
#else
extern graph G;
#endif
#endif // End of graph.h








