void augment(simpleQueue &Q)
{
  node *v;

  while (!Q.isEmpty()) {
    v = Q.pop();
    v->getCurrent()->use(); 
#ifdef STATS
    stats.pushCnt++;
#endif
  }

  G.getSource()->getCurrent()->use();
#ifdef STATS
    stats.pushCnt++;
#endif
}   

void  DFSUnitFlow(simpleQueue &Q, simpleQueue &S)
  // Q is used as DFS stack
  // S contains all visited nodes for reinitialization
{
  node *u, *v;
  arc  *a, *e, *stopA;
  bool pathFound;

  forallNodes(u,G) 
    {
      u->setReached(false);          // mark all nodes as unreached
      u->resetCurrent();
    }
  G.getSource()->setReached(true);
  forallOutArcs(a,G.getSource())
    {
      u = a->head();
      if (u->isReached()) continue;

      G.getSource()->setCurrent(a);

      if (u == G.getSink())
	{
	  augment(Q);
	  continue;
	}

      // DFS from u
      pathFound = false;
      Q.enqueue(u);
      S.enqueue(u);
      u->setReached(true);

      do
	{      
	  v = Q.pop();
	  e = v->getCurrent();

	  stopA = v->lastOutArc();
	  while (e <= stopA) // search forward
	    {
	      if (e->isUsed()) 
		{
		  e++; continue;
		}
	      u = e->head();
	      if (u->isReached())
		{
		  e++; continue;
		}

	      v->setCurrent(e);
	      Q.enqueue(v);

	      if (u == G.getSink()) 
		{
		  augment(Q);
	  	  pathFound = true;

		  break;
		}
	      else
		{
		  u->setReached(true);
		  S.enqueue(u);
		  v = u;
		  e = v->getCurrent();
		  stopA = v->lastOutArc();
		}
	    }
#ifdef STATS
	  stats.relabelCnt++;
#endif
	} while (!Q.isEmpty() && !pathFound);

      if (pathFound)
	{
	  // clean up
	  Q.reinitialize();
	  while (!S.isEmpty())
	    {
	      u = S.pop();
	      u->setReached(false);
	      u->resetCurrent();
	    }
	}
      else
	S.reinitialize(); // visited nodes stay reached
    }
}
