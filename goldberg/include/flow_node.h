// Hey EMACS! This is -*- c++ -*- 
// $Id: flow_node.h,v 1.12 1997/08/15 14:48:42 mslevine Exp $ 

#ifndef FLOW_NODE_H
#define FLOW_NODE_H

class graph;
#include <assert.h>
#include <constant.h>
#include <macros.h>
#include <node.h>

/***********************************************************************/
/*                                                                     */
/*         CLASS node_flow_data                                        */
/*                                                                     */
/***********************************************************************/

class node_flow_data : public node_basic_data
{
#define NOT_LABELED LARGE
  arc *current;
#if defined(BUCKET)
  node *bNext;     // next node in bucket
  node *bPrev;     // previous node in bucket 
#endif

public:

  operator long();
#ifndef BDFS
  void unlabel()                 {setDistance(NOT_LABELED);} 
  bool_t isUnlabeled()             {return (getDistance() == NOT_LABELED);}
#endif
  void setCurrent(arc *e)        {current = e;}
  void resetCurrent()            {current = firstOutArc();}
  arc  *getCurrent()             {return current;}
#if defined(BUCKET)
  node *getBNext()               {return bNext;}
  void setBNext(node *u)         {bNext = u;}
  node *getBPrev()               {return bPrev;}
  void setBPrev(node *u)         {bPrev = u;}
#endif
  void initialize()              {
				  resetCurrent();
				  node_basic_data::reinitialize();
#ifdef BUCKET
				  bNext = nil;
				  bPrev = nil;	
#endif
#ifndef BDFS
				  setDistance(0);
#else
                                  setReached(0);
#endif
  }

#ifdef NODE_CAPACITIES
  CapType getCapacity()           {return cap;}
  void setCapacity(CapType val)   {cap = val;}
  void decreaseCapacity(CapType val) {cap -= val;}
#endif // NODE_CAPACITIES


#ifdef PHASE2
  void setColor(long l)           {setDistance(l);}
#define WHITE 0
#define GREY  1
#define BLACK 2
  bool_t isWhite()                  {return getDistance()==WHITE;}
  bool_t isGrey()                   {return getDistance()==GREY;}
#endif // PHASE2

  CapType getOrigOutCapacity();
  CapType getOrigInCapacity();
};
#endif // End of flow_node.h



