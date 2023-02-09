/* GEOMETRIC SORTING */

void DiaSortItems(long *it, long N, double **d);
  /*
    Permutes "it[0..N-1]" so that "it[0]" and "it[N-1]"
    is a diameter, and the remaining items are sorted by 
    the relative distances to those two.
  */
  
void FarSortItems(long *it, long N, double **d);
  /*
    Sorts "it[0..N-1]" by repeatedly taking the element
    furthest from the path and inserting it in the best place.
  */
  
void CluSortItems(long *it, long N, double **d);
  /*
    Builds a dendrogram tree whose leaves are the items in "it[0..N-1]".
    Then swaps and/or reverses the children of every node so that 
    the jump between them is minimized. Finally enumerates the 
    leaves in this sorted tree and stores them back into "it". 
  */
  
void UniReSortItems(long *it, long N, double **d);
  /*
    Scans "it[0..N-1]" and tries to move each element
    so as to reduce the total path length. 
  */
  
void DelReSortItems(long *it, long N, double **d);
  /*
    Removes a big subset of the indices, and inserts them 
    back one by one. 
  */
  
long InsertItemInPath(long *it, long *Np, long ki, double **d);
  /*
    Inserts item "ki" in the best place of the path "it[0..*Np-1]".
    Assumes there is space in "it[]" for one more element.
    Increments *Np.
  */
  
