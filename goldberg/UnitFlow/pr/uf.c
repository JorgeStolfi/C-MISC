long relabel(node *u);
bool discharge(node *u);
void globalUpdate(simpleQueue &Q);

bucketArray ba;

void forwardUnitFlow(simpleQueue &Q)
{

  long relsSinceUpdate = 0; // used to determine when to perform global update
  node *u;

  long d=0, dT=0;
  arc *a;

  forallNodes(u,G) u->initialize();
  forallOutArcs(a,G.getSource()) d++;
  forallOutArcs(a,G.getSink()) dT++;
  d = min(d, dT);
  G.getSource()->setExcess(d);
  //  G.getSink()->setExcess(-d);
  G.getSink()->setExcess(-G.numArcs()-1); // so that the sink is never active
  globalUpdate(Q);

  while( (u=ba.getNextActiveNode()) != nil) {
    //    cout << u->index() << '\n';
    if (!discharge(u)) {             // u must be relabeled 
      relabel(u);
      relsSinceUpdate++;
#ifndef NO_BUCKET
      if (ba.checkGap())
	ba.removeGap();
#endif
      if ( ba.globalUpdateTime(relsSinceUpdate) ) {  
	globalUpdate(Q);
	relsSinceUpdate = 0;
      }
    }
#ifndef NO_BUCKET
    else {                              // u is inactive
      ba.makeInactive(u);
    }
#endif
  }
}

long relabel(node *u) 
{
  arc  *e;
  arc  *minA;  // the arc out of u with lowest head
  node *v;
  long minD;   // the height of the head 
  long d;

  assert(u->hasExcess());
  minD = LARGE;

#ifdef STATS
  stats.relabelCnt++;
#endif
  minA = u->lastOutArc();

  forallOutArcs(e,u) {
    if (e->isAvailable()) {
      d = e->head()->getDistance();
      if (d < minD) {
	minD = d;
	minA = e;
      }
    }
  }

  assert(minD+1 > u->getDistance());
  u->setDistance(minD+1);       // relabel u relative to minA
  if (minD+1 < ba.pathBound()) {
    u->setCurrent(minA);        // set the current arc pointer of u
    //cout << "Making Active "<< u->index() << "\n";
    ba.makeActive(u);           // make node u active again 
  }
  //cout << "Relabelling Node " << u->index() << " ";
  //cout << minD << "\n";
  return minD;
}
  
bool discharge(node *u)
{
  node *v;
  arc *e;
  long d;

  //cout << "Discharging Node " << u->index() << " ";
  assert(u->hasExcess());
  d = u->getDistance()-1;
  forallOutArcsAfter(e,u,u->getCurrent()) {
    v = e->head();
    if (v->getDistance() != d) continue;
    if (!e->isAvailable()) continue;
    //cout << v->index() << " ";
#ifdef STATS
    stats.pushCnt++;
#endif
    e->use();
    u->decrementExcess();
    if ( v->getExcess() == 0 ) {     // v is currently inactive
      ba.changeToActive(v);          // so v becomes active
      //cout << "Changing to Active "<< v->index() << "\n";
    }
    v->incrementExcess();
    if (!u->hasExcess()) {e++;break;}
  }
  u->setCurrent(e);
  //cout << e << "\n";
  return (u->hasExcess() ? false : true );
} 

void  globalUpdate(simpleQueue &Q)
{
  node *u, *v;
  arc  *e;
  DistType  currentDistance;             // the current depth in the breadth first search

#ifdef STATS
  stats.updateCnt++;
#endif

  Q.reinitialize();

  //  cout << "Performing Global Update\n";
  ba.reinitialize();                     // reinitialize the call bucketArray 
  forallNodes(u,G) u->increaseDistance(G.numNodes());          // mark all nodes as unlabeled
  G.getSink()->setDistance(0);           // set sink to distance 0
  Q.enqueue(G.getSink());                // add the sink to the queue

  while (!Q.isEmpty()) {      
    v = Q.dequeue();
    currentDistance = v->getDistance();  // distance of node being scanned

    forallOutArcs(e,v) {
      arc &ee = *e;
      u = ee.head();
      if (u->getDistance()-G.numNodes() < 0) continue;
      if (e->isAvailable()) continue;
      assert (currentDistance+1 >= u->getDistance()-G.numNodes());
      if (currentDistance+1 > u->getDistance()-G.numNodes()) u->resetCurrent();
      u->setDistance(currentDistance+1);
      //      cout << u->index() << ' ' << u->getDistance() << '\n';
      if ( u->hasExcess() ) {             // if u has excess make it active
	//cout << "Making Active "<< u->index() << "\n";
	ba.makeActive(u);                     
      }
      else {
	ba.makeInactive(u);
      }
      Q.enqueue(u);
    }
  }
}

#ifdef PHASE2
void phase2(simpleQueue &S, simpleQueue &topS)
{
  node *u, *v;
  node *root;
  arc *e;

  S.reinitialize();

  forallNodes(u,G) {
    u->setColor(WHITE);
    u->resetCurrent();
  }
  //  G.getSource()->setColor(BLACK);

  forallNodes(u,G) {
    if ((!u->isWhite()) || (!u->hasExcess()) || (u == G.getSink())) continue;
      
    u->setColor(GREY);
    root = u;
    do {
      while ((e = u->getCurrent()) <= u->lastOutArc()) {
	if (e->isAvailable() && !e->isForward() ) {
	  v = e->head();

	  if (v->isWhite()) {                    // start scanning v 
	    v->setColor(GREY);
	    S.enqueue(u);
	    u = v;
	    break;
	  }

	  else if (v->isGrey()) {                // cycle is found
	    do {                                 // remove cycle
	      u->getCurrent()->use();
	      u->setColor(WHITE);
	      u = S.pop();
	    } while (u != v);
	    u->setCurrent(u->getCurrent()+1);
	    /*
	    u = v;
	    do {                                 // remove cycle
	      v->getCurrent()->use();
	      if (v != u) v->setColor(WHITE);
	      v = v->getCurrent()->head();
	    } while (u != v); 
	    u->setCurrent(u->getCurrent()+1);
	    while (S.pop() != u) ;
	    */
	    break;
	    
	  }
	}
	e++;
	u->setCurrent(e);
      }
      if (u->getCurrent() > u->lastOutArc()) {         // scan of u complete
	u->setColor(BLACK);
	//	assert (u != G.getSource());
	topS.enqueue(u);
	if(u != root)
	  u = S.pop();
	else
	  break;
      }
    } while (true);
  }

  // return flow
  while (!topS.isEmpty()) {
    u = topS.pop();
    forallOutArcs(e,u) {
      if (!u->hasExcess() || u == G.getSource()) break;
      if (e->isAvailable() && !e->isForward()) {
	e->use();
	u->decrementExcess();
	e->head()->incrementExcess();
      }
    }
  }
}   
#endif	    
