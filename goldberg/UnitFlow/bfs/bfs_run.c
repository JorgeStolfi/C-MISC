// Hey emacs! This is -*- c++ -*-
// $Id: uf_run.cc,v 1.8 1997/08/15 17:51:09 pmartin Exp $

#define MAIN
#define UNIT_FLOW

#include <iostream.h>
#include <stdio.h>
#include <stdlib.h>

#include <param.h>
#include <statistics.h>
#include <uf.h>

#include <inline.h>

#include "timer.c"
#include "bfs.cc"
#include <../pr/uf_parser.h>
#include "base_parser.cc"
#include "parse_arguments.cc"

void BFSUnitFlow();
void parse();

main(int argc, char *argv[])
{
  parse_arguments(argc,argv,"go:u:");

  CapType flowValue;
  arc *e;

#ifndef CONCISE
  cout << "Begin Parsing...\n";
#endif
  parse();

#ifndef CONCISE
  cout << "Parsing has finished. Begin max-flow calculations...\n";
  cout << "\n";
  cout << "c Nodes: " << G.numNodes() << "\n";  
  cout << "c Arcs:  " << G.numArcs()/2 << "\n";
  cout << "\n";
#endif


  // need to do this alloc here
  simpleQueue Q(G.numNodes()); // BFS queue
  simpleQueue S(G.numNodes()); // stack containing nodes other than s with finite distance

  float t = timer();
  BFSUnitFlow(Q,S);
  t = timer() - t;

  flowValue = 0;
  forallOutArcs(e,G.getSource())
    if (e->isUsed())
      flowValue++;
  
#ifndef CONCISE
  cout << "c Running Time: " << t << '\n';
  cout << "c Flow: " << flowValue << '\n';
  cout << "c pushes: " << stats.pushCnt << '\n';
  cout << "c relabels: " << stats.relabelCnt << '\n';
#endif
#ifdef CONCISE
printf ("%8ld %8ld %8ld %10ld %10ld %11.2f\n",
         G.numNodes(),G.numArcs()/2,
	flowValue,
	stats.pushCnt,stats.relabelCnt,t);
#endif
}
