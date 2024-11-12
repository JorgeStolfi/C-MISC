
    /* Compute the coeff {C[N]} by the unit sum constraint: */
    double Csum = 0.0;
    for (int32_t k = k0; k < N; k++) { Csum += V[k]*p->C[k]; } 
    p->C[N] = (1.0 - Csum)/2;
    
    /* Set{C[0]} is fixed at zero: */
  }


/*  
      Let {R = N-k0} be the number of free coefficients, excluding {C[N]}.
      
      We can rewrite that constraint as 
      
        { 1 = 2*C[N] + SUM{ V[k]*C[k] : k \in k0..N-1 } }
        
      Hence we can compute {C[N]} from the other coefficients:
      
        { C[N] = 1 - (1/2) * SUM{ V[k]*C[k] : k \in k0..N-1 } }
      
      Note that this sets {C[N] = 1/2} if {N} is 1 and {noMiddle} is true.


  
      
      Let's rename {C[k0..N-1]} as {c[0..R-1]}; that is, {c[r]} is the same as {C[r+k0]},
      for all {r} in {0..R-1}.  We can then rewrite the approximant {F} as depending on 
      the "free" coefficients {c[0..R-1]} only:
      
        { F = G[N] + SUM{ c[r]*T[r] : r \in 0..R-1 } }
        
      where {T[r] = G[r+k0] - (V[r+k0]/V[N])*G[N]} for all {r} in {0..R-1}.
             
      The mismatch function {<F-H|F-H>} then too is a quadratic function 
      of the "free" coefficients {c[0..R-1]} only, namely
      
        { Q(c[0..R-1]) = 
            <G[N]-H|G[N]-H> + 
            2 * SUM{ c[r] * <G[N]-H|T[r]> : r \in 0..R-1 } + 
            SUM{ c[s]*c[t]*<T[s]|T[t]> : s,t \in 0..R-1 }
        }
       
      This function is minimized at a point {c[0..R-1]} when {dQ[r] = \partial_r Q = 0}
      for all {r} in {0..R-1}. This derivative is 
      
        { dQ[r] = 2*<G[N]-H|T[r]> + 2*SUM{ c[t]*<T[r]|T[t]> : t \in 0..R-1 }
        
      Therefore, we can determine the coefficients {C[k0..N-1]} by solving the first-degree
      equations {dQ[r] = 0} for all {r} in {0..R-1}.
      
      More precisely, eliminating the factor 2, the equations {dQ[r] = 0}
      become
      
        { SUM{ c[t]*<T[r]|T[t]> : t \in 0..R-1 } = - <G[N]-H|T[r]> }
      
      for each  {r} in {0..R-1}.  That is, we have the 
      linear equation system {U c = Y} where {U} is the  R×R}
      matrix with generic element { U[r][t] = <T[r]|T[t]>},
      {c} is considered a column vector, and {Y} is the column vector 
      with {Y[r] = - <G[N]-H|T[r]>}
        
      We then can solve the system to get {C[k0..N-1] = c[0..R-1]},
      compute {C[N]} from the unit-sum equation, and set {C[0] = 0}
      if {noMiddle} is true.

*/
