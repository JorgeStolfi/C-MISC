// Hey emacs! This is -*- c++ -*-
// $Id: bucket.h,v 1.10 1997/07/16 19:07:44 mslevine Exp $

#ifndef BUCKET_H
#define BUCKET_H

#ifdef QUEUE
#include <stack_queue.h>
#endif

#ifdef LIFO
#define Q_SELECT pop
#else
#define Q_SELECT dequeue
#endif

#if defined(BUCKET) && defined(NO_BUCKET)
#error you can only define one of BUCKET and NO_BUCKET
#endif
#if !(defined(BUCKET) || defined(NO_BUCKET))
#define NO_BUCKET
#endif

#if defined(QUEUE) && defined(NO_QUEUE)
#error you can only define one of QUEUE and NO_QUEUE
#endif
#if !(defined(QUEUE) || defined(NO_QUEUE))
#define NO_QUEUE
#endif


#ifdef BUCKET
#define INACTIVE_LIST 1
#define ACTIVE_LIST   2

//global vars
node *bucNextN, *bucPrevN, *bucU;

class bucketClass
{
  node *firstActive;
  node *firstInactive;
 
public:
  node *getFirstActive()              { return firstActive; }
  node *getFirstInactive()            { return firstInactive; }
  void setFirstActive(node *u)        { firstActive = u; }
  void setFirstInactive(node *u)      { firstInactive = u; }  
  bool_t inactiveListIsEmpty()          { return (firstInactive == nil); }
  bool_t activeListIsEmpty()            { return (firstActive == nil); }

  void addNode(node *u, int list)     { if (list == ACTIVE_LIST) addNode(u,firstActive);
                                        else addNode(u,firstInactive);}

  void addNode(node *u, node *&head)  { if (head != nil) head->setBPrev(u);
					u->setBNext(head);
					u->setBPrev(nil);
					head = u;
				      }

  void deleteNode(node *u, int list)  { if (list == ACTIVE_LIST) deleteNode(u,firstActive);
                                        else deleteNode(u,firstInactive);}

  void deleteNode(node *u, node *&head) {
					  bucNextN = u->getBNext();
					  if (head == u) {  // Deleting first element
					    head = bucNextN;
					    if (bucNextN != nil)
					      bucNextN->setBPrev(nil);
					  }
					  else {
					    bucPrevN = u->getBPrev();
					    bucPrevN->setBNext(bucNextN);
					    if (bucNextN != nil)
					      bucNextN->setBPrev(bucPrevN);
					  }
					}

  void removeLists()                   { removeList(firstActive);
					 setFirstActive(nil);
					 removeList(firstInactive);
					 setFirstInactive(nil);
				       }
  

  void removeList(node *&head)         {
					 for(bucU=head; bucU!=nil; bucU=bucU->getBNext()) {
					   bucU->setDistance(G.numNodes());
#ifdef STATS
					   stats.gNodeCnt++;
#endif
					 }
				       }

  void reinitialize()                  { firstActive=nil;
					 firstInactive=nil;
				       }

};
#endif  // BUCKET

class bucketArray
{
  DistType pthbnd;
  long n;
  long dMax;
#ifdef BUCKET
#ifndef QUEUE
  long aMax;
  long aMin;
#endif  // !QUEUE
  bucketClass *bucket;
  bucketClass *activeBucket;
#endif  // BUCKET
#ifdef QUEUE
  simpleQueue Q;
#endif  // QUEUE
  
public:
  bucketArray() {}
  bucketArray(graph &G) {initialize(G,G.numNodes());}
  bucketArray(graph &G, long l) {initialize(G,l);}
  void initialize(graph &G) {initialize(G,G.numNodes());}
  void initialize(graph &G, long l) {
#ifdef BUCKET
                              bucket = new bucketClass[l]; 
#endif // BUCKET
#ifdef QUEUE
			      Q.initialize(G.numNodes());
#endif // QUEUE
			      pthbnd=l; 
			      n=G.numNodes();
			      reinitialize(); }  

#ifdef BUCKET
  ~bucketArray()      { delete bucket; }
#endif // BUCKET

  globalUpdateTime(long relsSinceUpdate) { 
    return (relsSinceUpdate > param.globalUpdateFreq *(double)n );}

  bool_t checkGap()            { 
#ifdef BUCKET
                               return param.gap &&
				 (activeBucket->inactiveListIsEmpty() 
				  && activeBucket->activeListIsEmpty());
#else
			       return false;
#endif // BUCKET
			     }
  void reinitialize()        {
#ifdef BUCKET
                              long i; 
                              for(i=0; i<pathBound(); i++) 
				bucket[i].reinitialize();
#ifndef QUEUE
                              aMin=LARGE;aMax=dMax=0;
#endif // !QUEUE
#endif // BUCKET
#ifdef QUEUE
#ifndef NO_REBUILD
			      Q.reinitialize();
#endif
#endif
                             }

  void setPathBound(DistType l)  {pthbnd = l;}
  long pathBound()           {return pthbnd;}

  node *getNextActiveNode()  { node *u=NULL;
#if defined(BUCKET) && !defined(QUEUE)
			       if (aMax < aMin) return nil;
#endif
#if defined HI_SELECT
				 while ((u = bucket[aMax].getFirstActive()) == nil) {
				   aMax --;
				   if (aMax < aMin) return nil;
				 }
				 bucket[aMax].setFirstActive(u->getBNext());
				 activeBucket = bucket + aMax;
#elif defined LO_SELECT
				 while ((u = bucket[aMin].getFirstActive()) == nil) {
				   aMin ++;
				   if (aMax < aMin) return nil;
				 }
			       bucket[aMin].setFirstActive(u->getBNext());
			       activeBucket = bucket + aMin;
#elif defined HILO_SELECT
			       static last = 0;
			       last = !last;
			       if (last)
				 {
				   while ((u = bucket[aMax].getFirstActive()) == nil) {
				     aMax --;
				     if (aMax < aMin) return nil;
				   }
				   bucket[aMax].setFirstActive(u->getBNext());
				   activeBucket = bucket + aMax;
				 }
			       else
				 { 
				   while ((u = bucket[aMin].getFirstActive()) == nil) {
				     aMin ++;
				     if (aMax < aMin) return nil;
				   }
				   bucket[aMin].setFirstActive(u->getBNext());
				   activeBucket = bucket + aMin;
				 }
#elif defined QUEUE
			       do {
				 if (Q.isEmpty()) return nil;
				 u = Q.Q_SELECT();
			       } while(u->getDistance() >= pathBound());
#ifdef BUCKET
			       activeBucket = bucket + u->getDistance();
			       activeBucket->deleteNode(u,ACTIVE);
#endif // BUCKET
#endif // selection rule
			       return u;}

  void makeActive(node *u)   { 
#ifdef BUCKET
                               long d = u->getDistance();        // Depth of u
			       assert (d < pathBound());
			       bucket[d].addNode(u,ACTIVE_LIST);
#ifndef QUEUE
			       if (d < aMin) aMin=d;
			       if (d > aMax) aMax=d;
#endif // ! QUEUE
			       if (d > dMax) dMax=d;
#endif // BUCKET
#ifdef QUEUE		       
			       Q.enqueue(u);
#endif // QUEUE 
			     }
#ifdef BACKWARD_REBUILD
  void addFront(node *u) { Q.addFront(u);}
#endif

  void changeToActive(node *u)   { 
#ifdef BUCKET
                                   long d = u->getDistance();   // Depth of u
			           assert (d < pathBound());
				   bucket[d].deleteNode(u,INACTIVE_LIST); 
				   bucket[d].addNode(u,ACTIVE_LIST);
#ifndef QUEUE
				   if (d < aMin) aMin=d;
				   if (d > aMax) aMax=d;
#endif // ! QUEUE
#endif // BUCKET
#ifdef QUEUE		       
#ifdef PUSH_FRONT
				   Q.addFront(u);
#else
				   Q.enqueue(u);
#endif 
#endif // QUEUE
                                 }

  void move(node *u, DistType newD) {
#ifdef BUCKET
				   long d = u->getDistance();
				   if (d != newD)
				     {
				       int list;
				       if (u->hasExcess()) list = ACTIVE_LIST;
				       else list = INACTIVE_LIST;
				       bucket[d].deleteNode(u, list);
				       if (newD < pathBound())
					 {
					   bucket[newD].addNode(u, list);
#ifndef QUEUE
					   if (d < aMin) aMin=d;
					   if (d > aMax) aMax=d;
#endif // ! QUEUE
					   if (d > dMax) dMax=d;
					 }
				     }
#endif // BUCKET
				 }

  void makeInactive(node *u) {
#ifdef BUCKET
                               long d = u->getDistance();
			       assert (d < pathBound());
			       bucket[d].addNode(u,INACTIVE_LIST);
#endif
			     }


  void removeGap()           { 
#ifdef BUCKET			      
                               long r = activeBucket - bucket -1;
			       for(activeBucket++;activeBucket<=bucket+dMax;activeBucket++) {
				 activeBucket->removeLists();
			       }
			       dMax = r;
#ifndef QUEUE
			       aMax = r;
#endif // ! QUEUE
#endif // BUCKET
			     }
};

#endif // End of bucket.h
