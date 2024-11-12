void augment()
{
  node *v;
  arc *e;

  v = G.getSink(); 
  do 
    {
      e = v->getCurrent(); 
      e->unuse();
#ifdef STATS
      stats.pushCnt++;
#endif
      v = e->head();
    } 
  while (v != G.getSource());
}   

void  BFSUnitFlow(simpleQueue &Q, simpleQueue &S)
{
  node *u, *v;
  arc  *a, *e;
  bool foundPath;

  forallNodes(u,G) u->setReached(false); // mark all nodes as unlabeled
  G.getSource()->setReached(true);           // source is always reached
  forallOutArcs(a, G.getSource())
    {
      assert(!a->isUsed());
      u = a->head();
      if (u->isReached()) continue;

      // start BFS from u
      foundPath = false;
      u->setCurrent(a->reverseArc());
      if (u == G.getSink()) 
	{
	  foundPath = true;
	  augment();

	  continue;
	}

      u->setReached(true);
      Q.enqueue(u);
      S.enqueue(u);

      do 
	{      
	  v = Q.dequeue();
    
	  forallOutArcs(e,v) {
	    if (e->isUsed()) continue;
	    u = e->head();
	    if (u->isReached()) continue;
	    u->setCurrent(e->reverseArc());
	    assert(u->getCurrent()->head() == v);
	    if (u == G.getSink()) 
	      {
		foundPath = true;
		augment();

		break;
	      }
	    else
	      {
		u->setReached(true);
		Q.enqueue(u);
		S.enqueue(u);
	      }
	  }
#ifdef STATS
	  stats.relabelCnt++;
#endif
	} while (!Q.isEmpty() && !foundPath);

      // get ready for a new search
      Q.reinitialize();
      if (!foundPath)
	{
	  S.reinitialize(); // reached vertices remain marked
	}
      else
	{
	  // clean up
	  while (!S.isEmpty())
	    {
	      u = S.pop();
	      u->setReached(false);
	    }
	}
    }
}
