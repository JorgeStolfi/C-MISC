// -*- c++ -*-
// $Id: graph.inline.h,v 1.3 1997/06/20 13:40:25 mslevine Exp $
#ifndef GRAPH_INLINE_H
#define GRAPH_INLINE_H

#ifdef GRAPH_H
#include <assert.h>
/***********************************************************************/
/*                                                                     */
/*         inline functions for CLASS graph                            */
/*                                                                     */
/***********************************************************************/

inline void graph::makeNodes(long num, long alloc) {
  n_orig = n = num; n_alloc = alloc;
  node_array = new node[alloc + 1]; // need one node for sentinel
}

inline void graph::makeArcs(long num, long alloc) {
  m_orig= m = num; 
  m_alloc = alloc;
  arc_array = new arc[alloc + 1];  // need one arc for sentinel
}

inline long graph::name(node *v)               {return v-node_array+param.offset;}
inline long graph::index(node *v)              {return v-node_array;}
inline long graph::index(arc *e)               {return e-arc_array;}
inline node *graph::firstNode()                {return node_array;}
inline node *graph::lastNode()                 {return firstNode()+n-1;}
inline node *graph::succNode(node *v)          {return  v+1;}   
inline node *graph::precNode(node *v)          {return  v-1;} 
inline node *graph::node_i(long i)             {return node_array+i;} 

inline node *graph::newNode() {
  assert(n+1 < n_alloc); 
  return node_array+n++;
}

inline arc *graph::newArc() { 
  assert(m+1 < m_alloc);
  return arc_array+m++; 
}

inline arc *graph::firstArc()                  {return arc_array;}
inline arc *graph::lastArc()                   {return arc_array+m-1;}
inline arc *graph::succArc(arc *e)             {return e+1;}
inline arc *graph::precArc(arc *e)             {return e-1;}
inline arc *graph::arc_i(long i)               {return arc_array+i;} 

inline void graph::setSource(long i)           {src=node_i(i);}
inline void graph::setSource(node *v)          {src=v;}
inline node *graph::getSource()                {return src;}

inline void graph::setSink(long i)             {snk=node_i(i);}
inline void graph::setSink(node *v)            {snk=v;}
inline node *graph::getSink()                  {return snk;}

#endif // GRAPH_H

#endif // End of graph.inline.h
