// -*- c++ -*-
// $Id :$

#ifndef BFS_H
#define BFS_H

#include <assert.h>

typedef long  CapType;
class match_arc;
class match_node;
typedef match_node node;
typedef match_arc arc;

#include <constant.h>
#include <stack_queue.h>
#include <macros.h>
#include <node.h>
#include <graph.h>

#define FLOW_FORMULATION
#define UNMATCHED G.getSink()

typedef node Ynode;
typedef node Xnode;

class match_node
   : public node_basic_data
{
#define NOT_LABELED LARGE
#define LABELED 0
  node *mtch;
  node *current;

public:
  operator long();
  bool_t isMatched();
  node *getMatchNode()           {return mtch;}
  void match(node *v)            {mtch = v;}
  void unmatch();
  void setCurrent(node *v)       {current = v;}
  node  *getCurrent()            {return current;}
  void initialize(); 
};

class match_arc
{
//  static graph *G;
  node *hd;
public:
//  static void setGraph(graph *grph) {G=grph;}
  void setHead(node *v)             {hd=v;}
  node *head()                      {return hd;}
};

inline void match_node::initialize() {
  mtch = UNMATCHED;
  node_basic_data::reinitialize();
}

inline  bool_t match_node::isMatched()  {
return (mtch != UNMATCHED);
}

inline void match_node::unmatch()  {mtch = UNMATCHED;}

#define forallXNodes(u,G) \
    for( {arc *bfs_ee=(G.getSource())->firstOutArc(); \
         arc *bfs_stopA=G.getSource()->lastOutArc(); \
         u=bfs_ee->head();} \
	 bfs_ee <= bfs_stopA;\
	 u = (++bfs_ee)->head())

#define forallYNodes(u,G) \
    for( {arc *bfs_ee=(G.getSink())->firstOutArc(); \
         arc *bfs_stopA=G.getSink()->lastOutArc(); \
         u=bfs_ee->head();} \
	 bfs_ee <= bfs_stopA;\
	 u = (++bfs_ee)->head())

#endif // End of bim.h
