// -*- c++ -*-

#ifndef UF_H
#define UF_H

typedef long  CapType;
typedef long DistType;
class node_flow_data;
typedef node_flow_data node;
class unit_flow_arc;
typedef unit_flow_arc arc;

#include <assert.h>
#include <constant.h>
#include <stack_queue.h>
#include <macros.h>
#include <node.h>
#include <arc.h>
#include <flow_node.h>
#include <graph.h> 

class unit_flow_arc : public  arc_basic_data
{
  unit_flow_arc  *reverse;          // Pointer to reverse arc
  bool_t used;                        // true if already in use in this direction
  bool_t forward;                   // true is original input arc
#define REVERSE_ARCS                // Needed for correct parsing  
#define REVERSE_POINTER             // Needed for correct parsing  
public:
  bool_t isAvailable()                {return !used;}
  bool_t isUsed()                     {return used;}
  void setReverse(unit_flow_arc *e) {reverse = e;}
  arc *reverseArc()                 {return (unit_flow_arc *) reverse;}
  void use()                        {used=true; reverse->used=false;}
  void unuse()                      {used=false; reverse->used=true;}
  void setUsed(bool_t b)              {used = b;}      
  bool_t isForward()                  {return forward;}
  void setForward(bool_t b)           {forward = b;}
};

#endif \\  End of uf.h


