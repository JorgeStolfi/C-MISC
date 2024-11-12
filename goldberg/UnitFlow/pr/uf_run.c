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

#include <bucket.h>
#include <inline.h>

#include "timer.c"
#include "uf.cc"
#include <uf_parser.h>
#include "base_parser.cc"
#include "parse_arguments.cc"

void forwardUnitFlow(simpleQueue &Q);
void phase2();
void parse();

main(int argc, char *argv[])
{
  parse_arguments(argc,argv,"go:u:");

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
  simpleQueue Q(G.numNodes()); // for global update and also main stack for DFS

  ba.initialize(G);

  float t = timer();
  forwardUnitFlow(Q);
  float cut_time = timer() -t;
#ifdef PHASE2
  simpleQueue topS(G.numNodes()); // aux stack for DFS
  phase2(Q, topS);
#endif
  t = timer() - t;
#ifndef CONCISE
  cout << "c Running Time: " << t << '\n';
  //  cout << "c Flow: " << G.getSink()->getExcess() + min(G.getSink()->degree(),G.getSource()->degree()) << '\n';
  cout << "c Flow: " << G.getSink()->getExcess() +G.numArcs()+1 << '\n';
  cout << "c cut time: " << cut_time << '\n';
  cout << "c pushes: " << stats.pushCnt << '\n';
  cout << "c relabels: " << stats.relabelCnt << '\n';
  cout << "c updates: " << stats.updateCnt << '\n';
  cout << "c gap nodes: " << stats.gNodeCnt << '\n';
#endif
#ifdef CONCISE
arc *e;
printf ("%8ld %8ld %8ld %10ld %10ld %11.2f\n",
         G.numNodes(),G.numArcs()/2,
	//	min(G.getSink()->degree(),G.getSource()->degree())-
	//	G.getSource()->getExcess(), 
	G.getSink()->getExcess() + G.numArcs()+1,
	stats.pushCnt,stats.relabelCnt,t);
#endif
}
