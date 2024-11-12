// -*- c++ -*-

#ifndef UF_H
#define UF_H

#define WHITE 0
#define GREY 1
#define BLACK 2

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
  bool_t used;                        // true if already in use in this direction
  unit_flow_arc  *reverse;          // Pointer to reverse arc
#define REVERSE_ARCS                // Needed for correct parsing  
#define REVERSE_POINTER             // Needed for correct parsing  
  bool_t forward;
public:
  bool_t isUsed()                     {return used;}
  void setReverse(unit_flow_arc *e) {reverse = e;}
  arc *reverseArc()                 {return (arc *) reverse;}
  void use()                        {used=true; reverse->used=false;}
  void unuse()                      {used=false; reverse->used=true;}
  void setUsed(bool_t b)              {used = b;}      
  bool_t isForward()                  {return forward;}
  void setForward(bool_t b)           {forward = b;}
};

#endif \\  End of uf.h


