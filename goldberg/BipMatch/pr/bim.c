/* bim.c
 * Push-Relabel Bipartite Matching Algorithm
 *  
 * Requires classes:
 *
 * graph
 * 
 * Xnode *u
 *   XgetDistance()    - retutn the label of node v
 *   XsetDistance(D)   - set the label of node v to D    
 *   unlabel()         - mark node as unlabeled
 *   isUnlabeled()     - return true if v is unlabeled OR false otherwise
 *   initialize()      - reset node data
 *   getMatchNode()    - returns node matched to OR nil if unmatched
 *   match(Ynode *v)   - match nodes u and v. 
 *   unmatch(Ynode *v) - unmatch nodes u and v.
 *
 * Ynode *v
 *   YsetDistance(D)  - set the label of node v to D
 *   YgetDistance()   - return the label of node v
 *   isMatched()      - returns boolean indicating whether node is matched
 *   match(Xnode *u)  - match nodes u and v. 
 *   getMatchNode()    - returns node matched to OR nil if unmatched
 *   initialize()     - reset node data
 *
 * arc *e
 *   head()           - the head of arc e 
 *   isArtificial()   - return true if e is artificial false otherwise
 *                    - this is used to hide arcs from source and sink in max-flow
 *                      like formulations
 *  
 * bucketArray
 *   bucketArray(long l)     - create structure, l is the maximum posible path length. 	
 *   getNextActive()         - the next active node to be considered
 *   makeActive(node *v)     - make node v active
 *   changeToActive(node *v) - make node v, which is currently inactive, active 
 *   makeInActive(node *v)   - make node v inactive
 *   setPathBound(bound)     - set limit on path lengths
 *   pathBound()             - the current limit on path lengths
 *   checkGap()              - returns true if a gap has been found which should
 *                             be removed with removeGap
 *   removeGap()             - remove all nodes above the gap found 
 *   globalUpdateTime()      - return true if a global update should be performed
 *   reinitialize()          - reinitialize the class. Called initially and whenever
 *                             a global update is performed
 *      
 * Requires Macros: 
 *   forallXNodes(Xnode *u, graph *G)
 *   forallYNodes(Ynode *u, graph *G)
 *   forallOutArcs(arc *e, Xnode *u)              - all arcs out of X-node u
 *   forallOutArcs(arc *e, Ynode *v)              - all arcs out of Y-node v
 *   forallOutArcsAfter(arc *e, arc *f, Xnode *u) - all arcs out of X-node u occuring after 
 *                                                  arc f in the adjacency list
 *  ** All macros must assign nil to the loop variable upon normal loop termination **
 *
 */  

#include <math.h>
// Function Prototyptes
inline long relabel(Xnode *u);
void matchNode(Xnode *u);
void globalUpdate(simpleQueue &Q);
void greedyMatch();
bool augStage();

//  Xnode *u;     
//  Ynode *v, *x; 
//  arc  *e;

bucketArray ba;

void bim(simpleQueue &Q)
{
  long relsSinceUpdate = 0;             // used to determine when to perform global update
  long oldMVal = 0;
  Xnode *u;
  Ynode *v;
  forallYNodes(v,G) v->initialize();
  forallXNodes(u,G) {
    u->initialize();
    if (!param.greedy) ba.makeActive(u);     // initially all X-nodes are active
  }
if (param.greedy) greedyMatch();

#ifdef BFS
  if (augStage()) return;
#ifndef GREEDY
  globalUpdate(Q);
#endif
#endif

if (param.greedy) globalUpdate(Q);

#ifdef FLOW_FORMULATION
  G.getSink()->setDistance(-1);
#endif

  while( (u=ba.getNextActiveNode()) != nil) {

#if DEBUG_PRINT
    cout << "Active Node: " << u->index() << "\n";
#endif
    matchNode(u);                     
    if (!u->isMatched()) {              // u must be relabeled
      relabel(u);
      relsSinceUpdate++;
#ifndef NO_BUCKET
      if (ba.checkGap())
	ba.removeGap();
#endif
      if ( ba.globalUpdateTime(relsSinceUpdate) ) {
#ifdef AUG_STAGE
	  if (stats.mVal - oldMVal <= SWITCH_VAL) {
	    if (augStage()) return;
	  }
	  oldMVal = stats.mVal;
#endif
	globalUpdate (Q);
	relsSinceUpdate = 0;
      }
    }  
#ifndef NO_BUCKET                                 
    else {                              // u is matched and hence no longer active
      ba.makeInactive(u);
    }
#endif
  }
}  

void greedyMatch()
{
  Xnode *u;
  Ynode *v;
  arc *e;

  forallXNodes(u,G) {
    assert(!u->isMatched());
    forallOutArcs(e,u) {
      v = e->head();
#ifdef FLOW_FORMULATION
      if (v==G.getSource()) continue;
#endif 
      
      if (!v->isMatched()) {
	stats.mVal++;
	stats.greedCnt++;
	v->match(u);
	u->match(v);
	break;
      }
    }
    u->setCurrent(e);
  }
}

inline long relabel(Xnode *u)
{
  arc  *e;
  arc  *minA;  // the arc out of u with lowest head
  node *v;
  long minD;   // the height of the head 
  long d;

  minD = ba.pathBound();
#if DEBUG_PRINT
  cout << "Relabelling Node" << u->index() << " ";
#endif
  assert(!u->isMatched()); 

#ifdef STATS
  stats.relabelCnt++;
#endif

  minA = nil;

  forallOutArcs(e,u) {
    v = e->head();
#ifdef FLOW_FORMULATION
    if (v==G.getSource()) continue;
#endif 
    d = v->YgetDistance();
    if (d < minD) {
      minD = d;
      minA = e;
    }
  }

  u->XsetDistance(minD);        // relabel u relative to minA
  if (minD < ba.pathBound()) {
    u->setCurrent(minA);        // set the current arc pointer of u
    ba.makeActive(u);           // make node u active again  
  }
#if DEBUG_PRINT
  cout << minD << "\n";
#endif
  return minD;
}

void matchNode(Xnode *u)
{
  arc *e;
  Ynode *v;                      // v is a Y-node adjacent to u
  Xnode *w;                      // w is the X-node matched to v 

  assert(!u->isMatched());

  forallOutArcsAfter(e,u,u->getCurrent()) {
    v = e->head();
#ifdef FLOW_FORMULATION
    if (v == G.getSource()) continue;
#endif
    w = v->getMatchNode();
  
    if (w == UNMATCHED) {            // v is unmatched - so match u and v
#if DEBUG_PRINT
      cout << "Match Count Increased " << v->index() << "\n";
#endif

#ifdef STATS
      stats.mVal++;
      stats.pushCnt++;
#endif

      v->match(u);
      u->match(v);
      break;
    }

    else if (w->XgetDistance() < u->XgetDistance()) {  // w is lower then u

#if DEBUG_PRINT
      cout << "Double Push \n";
#endif

#ifdef STATS
      stats.pushCnt += 2;                // push from u to v AND from v to w          
#endif
      v->match(u);
      u->match(v);
      w->unmatch();
      ba.changeToActive(w);              // w is now unmatched to becomes active 
      break;
    }
    /*    if (u->isMatched())                 
      {
	//	v->setYDistance();
	break;                           // when u finds a match break 
      }*/
  }
  u->setCurrent(e);                      // set the current arc to the last arc examined 
}


void globalUpdate(simpleQueue &Q)
{
  Xnode *u;     
  Ynode *v, *x; 
  arc  *e;
  DistType  currentDistance;     // the current depth in the breadth first search

  //  cout << "Entering Global Update\n";
  //  simpleQueue <Ynode *> Q(G.numNodes()); // a queue containing Y-nodes to be scanned

#ifdef STATS
  stats.updateCnt++;
#endif

#if DEBUG_PRINT
  cout << "Performing Global Update\n";
#endif
  ba.reinitialize();                     // reinitialize the call bucketArray 
  forallXNodes(u,G) u->unlabel();        // mark all X-nodes as unlabeled

#ifdef FLOW_FORMULATION
  G.getSink()->setDistance(-1);
#endif

  forallYNodes(v,G) {
    if (!v->isMatched()) {               // add all unmatched Y-nodes to the queue
      //      v->YsetDistance(0);
      Q.enqueue(v);
    }
  }

  while (!Q.isEmpty()) {
    v = Q.dequeue();
    currentDistance = v->YgetDistance();  // distance of node being scanned
    forallOutArcs(e,v) { 
      u = e->head(); 
      // skip node if u has already been found or if e is a matching edge
      if ( u->isUnlabeled() && u != v->getMatchNode()) {
	u->XsetDistance(currentDistance);
	if ( !u->isMatched() ) {             // if u is unmatched make it active
#ifndef NO_REBUILD
#ifdef BACKWARD_REBUILD
	  ba.addFront(u);
#else
	  ba.makeActive(u);
#endif
#endif                     
	}
	else {                               // otherwise make is active  
#ifndef NO_BUCKET
	  ba.makeInactive(u);                // and add its match to the search tree
#endif
	  x = u->getMatchNode();
	  //	x->YsetDistance(currentDistance+1);
	  Q.enqueue(x);
	}
      }
    }
  }
  //  cout << "Exiting Global Update\n";
}

#ifdef AUGMENT
bool augStage()
{
  simpleQueue q(G.numNodes());
  //  simpleQueue st(G.numNodes());
  Xnode *u,*w;
  Ynode *v,*x;
  arc *e, *ee;
  bool doAugment;
  //  long size, tot_size=0;
  long iteration = 0;

  forallXNodes(u,G) u->setDistance(iteration);
#ifdef FLOW_FORMULATION
  G.getSink()->setDistance(LARGE);
  //  cout << "  " << G.getSink()->index() << " " << G.getSink()->getDistance() <<"\n";
#endif

  forallYNodes(v,G) {
    //cout << v->index() << "\n";
    if (v->isMatched()) continue;
    q.reinitialize();
    //    st.reinitialize();
    q.enqueue(v);
    doAugment = false;
    // size=0;
    iteration++;

    while((!doAugment) && (!q.isEmpty())) {
      v = q.dequeue();

      //cout << " " << v->index() << "\n";

      forallOutArcs(e,v) {
	u = e->head();
	if (u->getDistance()<iteration) {
	  //cout << "  " << u->index() << " " << u->getDistance() <<"\n";
	  u->setDistance(iteration);
	  //	  st.push(u);
	  // size++;
	  u->setPred(v);

	  /*
	  if (!u->isMatched()) {
	    doAugment=true;
	    break;
	  }*/

	  
	  //cout << "   " << u->getMatchNode()->index() << "\n";

	  x = u->getMatchNode();

	  forallOutArcsAfter(ee,x,x->getCurrent()) {
	    if (!ee->head()->isMatched()) {
	      u = ee->head();
	      u->setPred(x);
	      doAugment=true;
	      break;
	    }
	  }
	  x->setCurrent(ee);
	  if (doAugment) break;

	  q.enqueue(x);
	}
      }                     // node v is scanned 

      if (doAugment) {
	//cout << "Augmenting\n";
	stats.augCnt++;
	while(u!=UNMATCHED) {
	  //cout << " " << u->index() << "\n";
	  v = u->getPred();
	  //cout << "  " << v->index() << "\n";
	  w = v->getMatchNode();
	  //cout << "   " << w->index() << "\n";
	  v->match(u);
	  u->match(v);
	  u=w;
	}
	stats.mVal++;

	/*	while(!st.isEmpty()) {
	  u=st.pop();
	  u->unlabel();
	}
	*/

	/*
	tot_size +=size;
	if (tot_size > log((double) G.numNodes())*stats.mVal) {
	  cout << "Size: "<<size<< " (" << stats.mVal <<")\n";
	  return false;
	}
	*/
      }
    }
    if (!doAugment) {
      q.reuse();
      while (!q.isEmpty()) {
	v = q.dequeue();
	forallOutArcs(e,v) {
	  u = e->head();
	  if (u->getDistance()==iteration) u->setDistance(LARGE);
	}
      }
    }


    //    if (!doAugment) cout << "Delete Size: "<<size<<"\n";
  }
  return true;
}
#endif
