/* N-DIMENSIONAL GEOMETRY */

void NormalizeProbs(double *a, long N, long R);
  /*
    Scales the vector "a" so that it has unit "sum".

    Bombs out if "a"'s sum is zero or any "a[i]" is negative,
    reporting "R" as the line number.
  */

void EstimateProbs(double *a, long N, long R);
  /*
    Assumes the "a[i]" are integer event counts, and
    converts them to estimated probabilities by the
    formula "a[i] = (a[i]+1)/(sum + N)" where "sum" 
    is the sum of all "a[i]".

    Bombs out if any "a[i]" is negative,
    reporting "R" as the line number.
  */

double L2ProbDist(double *a, double *b, long N);
  /*
    Computes the Euclidean distance between vectors "a[0..N-1]" and
    b[0..N-1]", asssuming they are normalized to unit sum.  The
    distance is divided by sqrt(2) to give a result in "[0_1]".
  */
    
double EarthMoverDist(double *a, double *b, long N);
  /*
    Computes the "earth mover's" distance between vectors "a[0..N-1]"
    and b[0..N-1]", assuming they have been normalized to unit sum.
    The result is a number in "[0_1]".
  */
  
