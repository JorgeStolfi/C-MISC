// -*- c++ -*-
// $Id: flow_node.inline.h,v 1.3 1997/08/04 20:15:53 mslevine Exp $
#ifndef FLOW_NODE_INLINE_H
#define FLOW_NODE_INLINE_H
#ifdef FLOW_NODE_H

#ifdef FLOW_ARC_H // include only if class node_flow_arc is defined
inline CapType node_flow_data::getOrigOutCapacity()       
{
  CapType ans=0; arc *a;
  forallOutArcs(a, this)  ans += a->getOrigCapacity();
  return ans;
}

inline CapType node_flow_data::getOrigInCapacity()       
{
  CapType ans=0; arc *a=NULL;
  forallInArcs(a, this)  ans += a->getOrigCapacity();
  return ans;
}

#endif // FLOW_ARC_H 

#endif // FLOW_NODE_H 

#endif // End of flow_node.inline.h
