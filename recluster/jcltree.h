/* Last edited on 1997-11-24 00:51:47 by stolfi */
/* See the authorship and copyright notice at the end of this file. */

#ifndef jcltree_h
#define jcltree_h

#include <jclbasic.h>

/* ORDERED CLUSTER TREES:
  
  The procedures SortClusterTree and DumpClusterTree below take an
  integer N and two vectors "uch[0..N-2]" and "vch[0..N-2]", and interpret them 
  as an ordered binary trees with leaves "0..N-1" and internal nodes "N..2*N-2".
  
  A further signed parameter "R" specifies a subtree of that tree, and 
  a direction for its traversal, as follows: 
  
    If "R" is in the range "[0 .. N-1]", the subtree consists of the
    single leaf node number "R".

    If "R" is in the range "[N .. 2*N-2]", the subtree's root is the
    internal node number "R", whose left and right subtrees are
    "uch[R-N]" and "vch[R-N]", respectively.

    If "R < 0", the subtree in question is a left-to-right reversal of
    the subtree "-R".  In particular, if "R <= -N", the
    left and right subtrees are "-vch[-R-N]" and "-uch[-R-N]",
    respectively.

  The pointers "uch[]" and "vch[]" may be negative, in which case
  they too are interpreted according to the last rule above.
  
*/
  
void SortClusterTree(long *uch, long *vch, long R, long N, double **d, long *kap, long *kzp);
  /*
    Interprets "N", "uch", "vch" as an ordered binary tree on the
    items "0..N-1", as described above.  Traverses the subtree "R",
    swapping and/or reversing the two children of every internal node,
    if necessary, so as to minimize the length of the jump between
    them.  Also returns in "*kap" and "*kzp" the first and last leaves
    of the sorted subtree.
  */

void DumpClusterTree(long *uch, long *vch, long R, long N, long *it, long *ip);
  /*
    Interprets "N", "uch", "vch" as an ordered binary tree on the
    items "0..N-1", as described above.  Traverses the subtree "R"
    in its internal order, and stores its leaves in "it[*ip...]", 
    incrementing "*ip".
  */

#endif
