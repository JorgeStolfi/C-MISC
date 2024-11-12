// Hey emacs! This is -*- c++ -*-
// $Id: stack_queue.h,v 1.2 1997/06/12 21:52:23 mslevine Exp $

#ifndef STACK_QUEUE_H
#define STACK_QUEUE_H

long queueIndx;

class simpleQueue
{
  typedef node* data_t;
  data_t       *array;
  long       size;
  long       first_pos;
  long       last_pos;

  long  increment(long &l) {long tmp=l; l = (l+1) % size; return tmp;}
public:
  simpleQueue(long l)            {first_pos=last_pos=0;size=l+1; array = new data_t[l+1];}
  simpleQueue()                  {first_pos=last_pos=0;}
  void initialize(long l)             {size=l+1; array = new data_t[l+1];}
  ~simpleQueue()                 {delete array;}
  void enqueue(data_t x)         {array[increment(last_pos)]=x;}
  void addFront(data_t x)        {array[(first_pos = (size+first_pos-1) % size)]=x;}
  data_t dequeue()               {return array[increment(first_pos)];}
  data_t pop()                   {return array[(last_pos = (size+last_pos-1) % size)];}
  bool_t isEmpty()                 {return first_pos == last_pos;}
  void reinitialize()            {first_pos=last_pos=0;}
};

#ifdef OLD_CODE
template <class data_t> class simpleQueue
{
  data_t     *array;
  long       size;
  long       first_pos;
  long       last_pos;

  long  increment(long &l) {long tmp=l; l = l+1 % size; return tmp;}
public:
  simpleQueue(long l)            {first_pos=last_pos=0;size=l; array = new data_t[l];}
  ~simpleQueue()                 {delete array;}
  void enqueue(data_t x) {array[increment(last_pos)]=x;}
  data_t dequeue()        {if(isEmpty()) return nil; 
			   return array[increment(first_pos)];}
  isEmpty()               {return first_pos == last_pos;}
};
#endif

template <class data_t> class simpleStack
{
  data_t  *array;
  long    size;
  long    top_pos;

public:
  simpleStack(long l) {top_pos = 0; size=l; array = new data_t[l];}
  ~simpleStack()      {delete array;}
  void push(data_t x) {array[top_pos++] = x;}
  data_t pop()        {return array[--top_pos];}
  bool_t isEmpty()      {return top_pos == 0;}
};

#endif
