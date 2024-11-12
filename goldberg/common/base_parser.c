/********************************************************************/
/*                                                                  */
/*  parse (...) :                                                   */
/*       1. Reads the graph problem in extended DIMACS format.      */
/*       2. Prepares internal data representation                   */
/*                                                                  */
/*  Requires Routines:                                              */
/*       1. readProblemLine                                         */
/*       2. readNodeDescription                                     */
/*       3. readArcDescription                                      */
/*                                                                  */
/********************************************************************/
       
// Files to be included:

#include <stdlib.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream.h>
#include <assert.h>

//#ifdef NO_ERROR_HANDLING
//#define throw(string)  {cerr << string << "\n"; exit(1);}
//#define try
//#endif

inline void swap(arc *e1, arc *e2)
{
  arc temp = *e1;
  *e1 = *e2;
  *e2 = temp;
#ifdef REVERSE_POINTER
  if (e1->reverseArc() != e1)
    {
      e1->reverseArc()->setReverse(e1);
      e2->reverseArc()->setReverse(e2);
    }
  else 
    {
      e1->setReverse(e2);
      e2->setReverse(e1);
    }
#endif
}

void parse()
{
  long   *arc_tail=NULL;          // Internal Array: tails of the arcs
  long   *arc_first=NULL;         /* Internal Array for holding:
                                     - node degree
                                     - position of the first outgoing arc */
  long   tail;                    // Tail of current arc

  long   arc_num;
  long   arc_new_num;
//  DistType   max_cost = 0;

  arc    *arc_current=NULL;       // Pointer to the arc structure
  arc    *arc_new;
  node   *v; 

  long    no_lines=0;             // # of current input line
  long    no_plines=0;            // # of problem-lines
  long    no_nlines=0;            // # of node(source)-lines
  long    no_alines=0;            // # of arc-lines

  char    input_line[MAXLINE];    // For reading input line

/* The main loop:
        -  reads the line of the input,
        -  analises its type,
        -  checks correctness of parameters,
        -  puts data to the arrays,
        -  does service functions
*/
  try{
    while ( gets ( input_line ) != NULL ) {
      no_lines ++;

      switch (input_line[0]) 
	{
	case 'c':                       // Skip lines with comments
	case '\n':                      // Skip empty lines   
	case '\0':                      // Skip empty lines at the end of file
	case 't':                       // Name of the problem
	  break;
	    
	case 'p':                       // Problem description     
	  if ( no_plines > 0 ) throw("more than one problem line"); 
	  no_plines = 1;
   
	  readProblemLine(G,input_line);

	  arc_tail = new long[G.numArcs()]; 
	  arc_first= new long[G.numNodes()+1];
	  for (long i=0; i<G.numNodes()+1; i++) arc_first[i]=0;

	  if ( arc_first == NULL || arc_tail == NULL )
	    throw("can't obtain enough memory to solve this problem");
		     
	  arc_current = G.firstArc();	    // Setting pointer to the current arc
	  break;
	
	case 'n':		               // Node description
	  if ( no_plines == 0 ) throw("problem description must come first");
	  if ( no_nlines > MAX_NODE_LINES) throw("too many nodes in the input");

	  no_nlines++;   
	  readNodeDescription(G,input_line);
	  break;

	
	case 'a':                    // Arc Description
	  if ( no_nlines < MIN_NODE_LINES )  throw("too few nodes lines"); 	    
	  if ( no_alines >= G.numArcs() )  throw("too many arcs input");
		
	  tail = readArcDescription(G,input_line,arc_current);
	  arc_first[tail + 1] ++; /* no of arcs outgoing from tail
                                           is stored in arc_first[tail+1] */

	  // Storing Information About The Arc
	  arc_tail[no_alines] = tail;

#ifdef REVERSE_ARCS  //Add reverse arc as well
	  tail = G.index(arc_current->head());
	  arc_current++;
	  no_alines++;
	  arc_first[tail + 1] ++;
	  arc_tail[no_alines] = tail;
#endif
	  no_alines ++;
	  arc_current ++;
	  break;

	default:
	  throw("unknown line type in the input");
	  break;
	  
	} // End of switch
    }     // End of input loop

/* ----- all is red  or  error while reading ----- */ 

    if ( feof (stdin) == 0 ) throw("reading error");

    if ( no_lines == 0 ) throw("input file is empty");

//if ( no_alines < G.numArcs() ) // Not enough arcs
//  parserError(EN19,no_lines);

#ifdef ARTIFICIAL_SOURCE
    //    G.setSource(SOURCE);

    forallNodes(v,G)  {
      if (!connectedToSource(v)) continue;

      createSourceArc(arc_current,v);  // Create reverse arc at same time
      arc_current++;

      tail =  G.getSource()->index();
      arc_first[tail + 1] ++;
      arc_tail[no_alines] = tail;
      no_alines++;
#ifdef REVERSE_ARCS      
      arc_current++;
      tail =  G.index(v);
      arc_first[tail + 1] ++;
      arc_tail[no_alines] = tail;
      no_alines++;
#endif // REVERSE_ARCS

    }
#endif // ARTIFICIAL_SOURCE

#ifdef ARTIFICIAL_SINK  
    //    G.setSink(SINK);

    forallNodes(v,G)  {
      if (!connectedToSink(v)) continue;

      createSinkArc(arc_current,v);  // Create reverse arc at same time
      arc_current++;

      tail =  G.index(v);
      arc_first[tail + 1] ++;
      arc_tail[no_alines] = tail;
      no_alines++;
#ifdef REVERSE_ARCS      
      arc_current++;
      tail =  G.getSink()->index();
      arc_first[tail + 1] ++;
      arc_tail[no_alines] = tail;
      no_alines++;
#endif // REVERSE_ARCS

    }
#endif // ARTIFICIAL_SINK

    forallNodes(v,G) {
      long indx;
      indx = G.index(v);
      if (v !=  G.firstNode())   arc_first[indx] += arc_first[indx-1];
      v->initAdjList( G.arc_i(arc_first[indx]) );
    }
    // init. sentinel
    G.node_i(G.numNodes())->initAdjList(G.arc_i(G.numArcs()));

    /*
    // setup first arcs, including sentinel node
    long indx, sum;
    sum = 0;
    for (indx = 0; indx <= G.numNodes(); indx++)
      {
	sum+=arc_first[indx];
	v = G.node_i(indx);
	v->initAdjList(G.arc_i(sum));
      }
      */
/********** ordering arcs - linear time algorithm ***********/

/* before below loop arc_first[i+1] is the number of arcs outgoing from i;
   after this loop arc_first[i] is the position of the first 
   outgoing from node i arcs after they would be ordered;
   this value is transformed to pointer and written to node.first[i]
   */


    forallNodes(v,G) {
      long indx;
      indx = G.index(v);
      // AVG      last =  G.index(v->lastOutArc());
                             /* arcs outgoing from i must be cited    
                              from position arc_first[i] to the position
                              equal to initial value of arc_first[i+1]-1  */
      arc *e;

      e = G.arc_i(arc_first[indx]);
      forallRemainingOutArcs(e,v) { 
	arc_num = G.index(e);
	tail = arc_tail[arc_num];
          /* the arc no  arc_num  is not in place because arc cited here
             must go out from indx;
             we'll put it to its place and continue this process
             until an arc in this position would go out from indx */
	while ( tail != indx ) {
	  arc_new_num  = arc_first[tail];
	  arc_new      = G.arc_i(arc_new_num);
	  swap(arc_new,e);
	  arc_tail[arc_num] = arc_tail[arc_new_num];

	  // We Increase arc_first[tail] but Label Previous Position

	  arc_tail[arc_new_num] = tail;
	  arc_first[tail] ++ ;
	  
	  tail = arc_tail[arc_num];
	}
      }
      /* all arcs outgoing from  v  are in place */
    }       

// -----------------------  Arcs Are Ordered  ------------------------- */

// Assigning Output Values

#ifdef SHORTEST_PATH
if ( G.getSource()->adjListIsEmpty() ) throw("no arc out of source");
#endif
   }
#ifndef NO_ERROR_HANDLING
catch(char *str) {
  cerr << input_line;
  cerr << "\nLine " << no_lines << " of input - " << str << "\n";
  abort();
}
#endif

// Free Internal Memory 
  delete arc_first ; 
  delete arc_tail ;

/* ---------------------------------- */

}

/* --------------------   end of parser  -------------------*/


