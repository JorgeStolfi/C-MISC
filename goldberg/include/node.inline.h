// -*- c++ -*-
// $Id: node.inline.h,v 1.1 1997/06/19 19:10:03 pmartin Exp $
#ifndef NODE_INLINE_H
#define NODE_INLINE_H
#ifdef NODE_H

inline node::operator long()
{
  return (this - G.firstNode());
}

/***********************************************************************/
/*                                                                     */
/*         FUNCTIONS FOR CLASS node_basic_data                         */
/*                                                                     */
/***********************************************************************/
// inline static void node_basic_data::setGraph(graph *grph) {G=grph;}


//******************

inline arc *node_basic_data::firstOutArc() 
{return adj_list;}

//******************
inline arc *node_basic_data::lastOutArc()  
{return ((((node *) this+1)->firstOutArc())-1);}

//******************

inline arc *node_basic_data::succOutArc(arc *e) 
{return e+1;} 

//******************

inline arc *node_basic_data::prevOutArc(arc *e) 
{return e-1;}

//******************

inline bool_t node_basic_data::contains(arc *e)
{return ( (e >= firstOutArc()) && (e <= lastOutArc()) );} 

//******************

inline void node_basic_data::reinitialize()
{ 
#ifndef BDFS
  setDistance(VERY_FAR);
#else
  setReached(0);
#endif

#ifdef EXCESSES
  excess = 0;
#endif
}

inline long node_basic_data::name()
{
  return G.name((node *) this);
}

inline long node_basic_data::index()
{
  return G.index((node *) this);
}

inline bool_t node_basic_data::isSource()
{
  return (this == G.getSource());
}

inline bool_t node_basic_data::isSink()
{
  return (this == G.getSink());
}

inline bool_t node_basic_data::isTerminus()
{
  return (isSink() || isSource());
}

#endif // NODE_H
#endif // End of node.inline.h
