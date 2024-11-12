// Hey emacs! This is -*- c++ -*-
// $Id: node.h,v 1.9 1997/07/03 14:02:33 mslevine Exp $

// BDFS stands for breadth/depth first search

#ifndef NODE_H
#define NODE_H 

#if defined(EXCESSES) && defined(NO_EXCESSES)
#error you can only define one of EXCESSES and NO_EXCESSES
#endif

#if !(defined(NO_EXCESSES) || defined(EXCESSES))
#define NO_EXCESSES
#endif



/***********************************************************************/
/*                                                                     */
/*         CLASS node_basic_data                                       */
/*                                                                     */
/***********************************************************************/

class node_basic_data
{
  arc      *adj_list;      // adjacency list - pointer into array $arcs$

#ifdef EXCESSES
  CapType excess;     // current excess at node
#endif

protected:
#ifndef BDFS
  DistType distance;
#else
  bool_t reached;
#endif
//  static graph *G;
  
public:
//  static void setGraph(graph *grph) {G=grph;}
  inline void reinitialize();

  inline arc *firstOutArc();
  inline arc *lastOutArc(); 
  inline arc *succOutArc(arc *e);
  inline arc *prevOutArc(arc *e);
  inline bool_t contains(arc *e);

  void initAdjList(arc *e) {adj_list = e;}
#ifndef BDFS
  DistType getDistance()        {return distance;}
  void setDistance(DistType l)  {distance = l;}
  void changeDistance(DistType l) {distance += l;}
  void increaseDistance(DistType l) {distance += l;}
#else
  bool_t isReached() {return reached;}
  void setReached(bool_t i) {reached = i;}
#endif

  long index();
  long name();

  bool_t isTerminus();
  bool_t isSink();
  bool_t isSource();

#ifdef EXCESSES
  bool_t hasExcess()               { return ((excess >0) ? true :false);}
  void increaseExcess(CapType x) {excess += x;}
  void incrementExcess()         {excess++;}
  void decreaseExcess(CapType x) {excess -= x;}
  void decrementExcess()         {excess--;}
  void setExcess(CapType x)      {excess=x;}
  long getExcess()               {return excess;}  
#else
  bool_t hasExcess()               {return false;}
  long getExcess()               {return 0;}
#endif

};

#endif
