// -*- c++ -*-
// $Id: inline.h,v 1.4 1997/06/19 22:01:58 pmartin Exp $
#ifndef INLINE_H
#define INLINE_H
#else
#warning File inline.h include multiple times
#endif

#ifdef ARC_H
#include <arc.inline.h>
#endif // ARC_H

#ifdef FLOW_ARC_H
#include <flow_arc.inline.h>
#endif 

#ifdef FLOW_NODE_H
#include <flow_node.inline.h>
#endif 

#ifdef NODE_H
#include <node.inline.h>
#endif

#ifdef GRAPH_H
#include <graph.inline.h>
#endif

#ifdef TREE_H
#include <tree.inline.h>
#endif
