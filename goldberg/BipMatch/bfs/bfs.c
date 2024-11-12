#include <math.h>

void greedyMatch()
{
  Xnode *u;
  Ynode *v;
  arc *e;

  forallXNodes(u,G) {
    assert(!u->isMatched());
    forallOutArcs(e,u) {
      v = e->head();
      
      if (!v->isMatched()) {
	stats.mVal++;
	stats.greedCnt++;
	v->match(u);
	u->match(v);
	break;
      }
    }
  }
}

void augment(Ynode *last)
{
  Ynode *v, *v1;
  Xnode *w;
  arc *e;

  v = last;
  do 
    {
      w = v->getCurrent(); 
#ifdef STATS
      stats.pushCnt++;
#endif
      v->match(w);
      if (w->isMatched())
	{
	  v1 = w->getMatchNode();
	  assert(v != v1);
	  w->match(v);
	  v = v1;
	}
      else
	{
	  w->match(v);
	  break;
	}
    } 
  while (1);
#ifdef STATS
  stats.mVal++;
#endif
}   

void bfsBIM(simpleQueue &Q, simpleQueue &S)
{
  Xnode *u, *w;
  Ynode *v;
  arc *e;
  bool foundPath;

  forallYNodes(v,G) v->initialize();
  forallXNodes(u,G) 
    {
      u->initialize();
      u->setReached(false);
    }      

  if (param.greedy) greedyMatch();

  forallXNodes(u,G) 
    {
      if (u->isMatched() || u->isReached()) continue;

      // start BFS from u
      foundPath = false;
      u->setReached(true);
      Q.enqueue(u);
      S.enqueue(u);
      
      do
	{
	  u = Q.dequeue();

	  forallOutArcs(e,u)
	    {
	      v = e->head();
	      if (v == G.getSource()) continue;
	      if (v->isMatched())
		{
		  w = v->getMatchNode();
		  if (w->isReached()) continue; // includes w == u
		  assert(w != u);
		  v->setCurrent(u);
		  w->setReached(true);
		  Q.enqueue(w);
		  S.enqueue(w);
		}
	      else
		{
		  foundPath = true;
		  v->setCurrent(u);
		  augment(v);
		
		  break;
		}
	    }
#ifdef STATS
	  stats.relabelCnt++;
#endif
	} while (!foundPath && !Q.isEmpty());

      // get ready to start new search
      Q.reinitialize();
      if (!foundPath)
	S.reinitialize(); //reached vertices remaine marked
      else
	{
	  // clean up
	  while (!S.isEmpty())
	    {
	      u = S.pop();
	      assert(u->isMatched());
	      u->setReached(false);
	    }
	}
    }
}

