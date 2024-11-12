#ifndef MACROS_H
#define MACROS_H


#define forallArcs(e,G) for(e=(G).firstArc(); e <= G.lastArc(); e++)

#define forallNodes(v,G) for(v=(G).firstNode(); v <= G.lastNode(); v++)

#define forallRemainingOutArcs(e,v) \
          for(arc *macStopA=v->lastOutArc(); e<=macStopA; e++)

#define forallOutArcsAfter(e,v,ee) for(e=ee; e <= v->lastOutArc(); e++)

//#define forallOutArcs(e,v) for(e=(v)->firstOutArc(); e <= v->lastOutArc(); e++)
#define forallOutArcs(e,v) \
  for({arc *macStopA=v->lastOutArc(); e=(v)->firstOutArc();}\
  e <= macStopA; e++)

#define forallInArcs(e,v) for(arc *__re=(v)->firstOutArc();  \
			      (__re <= v->lastOutArc()) &&  \
			      (e = __re->reverseArc()); \
			      __re++)


#define forallIncidentArcs(e,v) \
        for(e=(v)->firstOutArc(),int _dir=0; _dir==1 || e <= v->lastOutArc() ; \\
	      e = (dir == 0) ? e->reverseArc() : e->reverseArc()+1, dir = (dir+1)%2)

#define forallNodesDeeperInTree(u,v) \
  for(u = v->getTreeNext(); (u->getTreeDepth() > v->getTreeDepth()); u=u->getTreeNext())


#define forallArcsOfCycle(e,v) \
        for({node *_u; \
	    _u = nil; e=v->getParentArc();} _u != v; \
	    _u = e->head()->getParent(), e = ((arc *)e->reverse)->head()->getParentArc())


#define forallArcsInList(e,L) for({int i=0; e=L.item(i);} e != nil; i++, e=L.item(i))

template <class X>
inline X max(X a, X b)
{ return (a > b) ? a : b;}

template <class X>
inline X min(X a, X b)
{ return (a < b) ? a : b;}

template <class X>
inline X abs(X a)
{ return (a < 0) ? -a : a;}

#endif






