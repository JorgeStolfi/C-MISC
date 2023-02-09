void DiaSortItems(long *it, long N, double **d)
  {
    /* Sorting keys: */
    double *rd;
    long ka, kb;
    if (N <= 2) return;
    /* Find a diametral pair ka,kb of items, put it in it[0], it[N-1]: */
    { long i, j, iMax, jMax;
      double dMax = -1;
      for (i=1; i<N; i++)
        for (j=0; j<i; j++)
          { double dd = d[it[i]][it[j]]; 
            if (dd > dMax) { dMax = dd; iMax = i; jMax = j; }
          }
      ka = it[iMax]; it[iMax] = it[0];   it[0] = ka;
      kb = it[jMax]; it[jMax] = it[N-1]; it[N-1] = kb;
    }
    /* Define sort keys based on distances to ka and kb: */
    { long i;
      double eps = 0.000001 * d[ka][kb];
      if (eps == 0) return;
      rd = (double *)(Alloc(N*sizeof(double)));
      rd[0] = -1; rd[N-1] = +1;
      for (i=1; i<N-1; i++)
        { long ki = it[i];
          double da = d[ki][ka]; 
          double db = d[ki][kb];
          rd[i] = (da-db)/(da+db+eps);
          if ((rd[i] <= -1) || (rd[i] >= +1)) { Error("program error - rd[i] out of range"); }
        }
    }
    /* Sort remaining items "it[1..N-2]" by ratio of distances to those two: */
    { long i;
      for (i=1; i<N-1; i++)
        { long ki = it[i];
          double ri = rd[i];
          long j = i;
          while (rd[j-1] > ri) { rd[j] = rd[j-1]; it[j] = it[j-1]; j--; }
          if (i != j) { rd[j] = ri; it[j] = ki; }
        }
    }
    if (rd != NULL) free(rd);
  }

void FarSortItems(long *it, long N, double **d)
  {
    /* Distances from path: */
    long ka, kb;
    if (N <= 2) return;
    /* Find a diametral pair of items ka,kb, put them in it[0], it[1]: */
    { long i, j, iMax, jMax;
      double dMax = -1;
      for (i=1; i<N; i++)
        for (j=0; j<i; j++)
          { double dd = d[it[i]][it[j]]; 
            if (dd > dMax) { dMax = dd; iMax = i; jMax = j; }
          }
      ka = it[iMax]; it[iMax] = it[0]; it[0] = ka;
      kb = it[jMax]; it[jMax] = it[1]; it[1] = kb;
    }
    /* Repeately look for the farthest item, insert it in path: */
    { long i, j;
      double *dd = (double *)(Alloc(N*sizeof(double)));
      /* Initialize the distance-from-path and best-place tables: */
      /* If current path if "it[0..i-1]", then */
      /*   The candidates for insertion are "it[i..N-1]" */
      /*   Element "dd[j]" is distance of "it[j]" to the path *vertices* */
      dd[0] = 0;
      dd[1] = 0;
      for (j=2; j<N; j++) 
        { long kj = it[j];
          double da = d[kj][ka]; 
          double db = d[kj][kb];
          dd[j] = (da < db ? da : db);
        }
      /* Now expand that path: */
      i = 2;
      while(i < N)
        { long j, jMax, ki, where;
          double ddMax;
          /* Look for item "it[jMax]" with largest "dd[jMax]": */
          ddMax = -1;
          for (j=i; j<N; j++)
            { double ddj = dd[j]; 
              if (ddj > ddMax) { ddMax = ddj; jMax = j; }
            }
          /* Save that item in "ki" and move the current "it[i]" there: */
          /* We don't need "dd[jMax]" anymore: */
          ki = it[jMax]; 
          it[jMax] = it[i]; dd[jMax] = dd[i]; 
          dd[i] = 0;
          where = InsertItemInPath(it, &i, ki, d);
          for (j=i; j<N; j++)
            { double dij = d[ki][it[j]];
              if (dij < dd[j]) { dd[j] = dij; }
            } 
        }
      if (dd != NULL) free(dd);
    }
  }

void CluSortItems(long *it, long N, double **d)
  {
    /* The first pass builds a tree. */
    /* The nodes of the tree are numbered in the range "[0..2*N-2]". */
    /* If "i" is in "[0..N-1]", node "i" is a leaf, which is item "i". */
    /* If "i" is in "[N..2*N-2]", node "i" is an internal node, whose children */
    /* are nodes "uch[i-N]" and "vch[i-N]". */
    long M = 2*N-1;
    long *uch = (long *)Alloc((M-N)*sizeof(long));
    long *vch = (long *)Alloc((M-N)*sizeof(long));
    
    /* During construction we have K forests whose roots nodes are "root[0..K-1]": */
    long K;
    long *root = (long *)Alloc(N*sizeof(long));
    
    /* The distance between forest "f" and forest "g" is stored in "e[i][j]": */
    double **e = (double **)Alloc(N*sizeof(double*));
    
    /* Initialize root set: */
    { long i; for (i=0; i<N; i++) { root[i] = it[i]; } }
    K = N;
    
    /* Intialize cluster distance matrix: */
    { long i, j;
      for (i=0; i<N; i++) 
        { e[i] = (double*)Alloc(N*sizeof(double));
          for (j=0; j<N; j++) { e[i][j] = d[i][j]; }
        }
    }
    
    /* Repeatedly merge the two nearest clusters: */
    while (K > 1)
      { /* The internal nodes created so far are [N..M-K] */
        /* Find the two nearest clusters */
        long iMin, jMin;
        { double dMin = 999;
          long i, j;
          for (i=0; i<K; i++)
            for (j=i+1; j<K; j++)
              { double dij = e[i][j];
                if (dij < dMin) { dMin = dij; iMin = i; jMin = j; }
              }
        }
        /* Merge them into a new cluster: */
        { long r = M-K+1; long k;
          uch[r-N] = root[iMin]; 
          vch[r-N] = root[jMin];
          /* Store the new root in root[iMin], */
          /* and compute its distances to other subtrees: */
          root[iMin] = r;
          for (k=0; k<K; k++)
            { if ((k != iMin) && (k != jMin))
                { double dik = e[iMin][k]; double djk = e[jMin][k];
                  /* double drk = (dik + djk + (dik < djk ? dik : djk)) / 3; */
                  double drk = 2/(1/dik + 1/djk);
                  e[iMin][k] = drk; e[k][iMin] = drk;

                }
             }
          /* Eliminate root[jMin]: */
          root[jMin] = root[K-1];
          K--;
          for (k=0; k<K; k++)
            { if (k != jMin)
                { double duk = e[K][k]; e[jMin][k] = duk; e[k][jMin] = duk; }
            }

        }
      }

    /* Second pass: sort the subtrees by proximity */
    { long ka, kz; SortClusterTree(uch, vch, +(M-1), N, d, &ka, &kz); }
    
    /* Third pass: dump it out: */
    { long i = 0;
      DumpClusterTree(uch, vch, +(M-1), N, it, &i);
      if (i != N) { Error("program error: DumpClusterTree"); }
    }
    
    free(e); free(uch); free(vch);
  }

void DelReSortItems(long *it, long N, double **d)
  {
    long i;
    if (N <= 2) return;
    /* Split path into two subsets "A = it[0..i-1]" and "B = it[i..N-1]": */
    srand(62712357);
    { long j;
      bool_t keep = TRUE;
      double db = d[it[0]][it[1]];
      double da = db;
      i = 0;
      for(j=1; j<N; j++)
        { /* The B items are "it[i..j-1]" */
          /* Break preferentially at long edges: */
          if (j == N-1)
            { keep = TRUE; }
          else 
            { double dc = ((j > N-2) ? db : d[it[j]][it[j+1]]);
              double er = db/(da + db + dc + 1.0e-100);
              double prob = 0.125 + 0.375 * er;
              if (rand() < (int)(prob*32767)) 
                { keep = !keep; 
                  /* fprintf(stderr, "keep=%d at j=%ld\n", keep, j); */
                }
              da = db; db = dc;
            }
          if (keep) 
            { /* keep it in the path */
              long k = it[i]; it[i] = it[j]; it[j] = k; i++;
            }
        }
    }
    /* Insert the B items into the A path: */
    while (i<N)
      { long ki = it[i]; long where;
        where = InsertItemInPath(it, &i, ki, d);
      }
  }
  
long InsertItemInPath(long *it, long *Np, long ki, double **d)
  {
    long p, pMin;
    long N = (*Np);
    double ddMin;
    /* Find the best edge "(ka,kb) = (it[pMin-1],it[pMin])" where to insert "ki": */
    ddMin = 999;
    for (p=0; p<=N; p++)
      { long ka = (p > 0 ? it[p-1] : -1);
        long kb = (p < N ? it[p] : -1);
        double da = (ka > 0 ? d[ki][ka] : 0); 
        double db = (kb > 0 ? d[ki][kb] : 0);
        double de = ((ka > 0) && (kb > 0) ? d[ka][kb] : 0);
        double ddp = da + db - de;
        if (ddp < ddMin) { ddMin = ddp; pMin = p; }
      }
    /* Open up that edge: */
    p = N;
    while (p > pMin) { it[p] = it[p-1]; p--; }
    it[pMin] = ki;
    N++;
    /* Update path length: */
    (*Np) = N;
    return (pMin);
  }
  
void UniReSortItems(long *it, long N, double **d)
  {
    long i;
    if (N <= 2) return;
    /* Take each item and see if it would fit better elsewhere */
    i = 0;
    while (i<N)
      { long ki = it[i];
        long jMin;
        double ci, cjMin; 
        /* Compute "ci" = gain in removing node "it[i]" from the path: */
        if (i==0) 
          { ci = d[ki][it[1]]; }
        else if (i==N-1)
          { ci = d[ki][it[N-2]]; }
        else
          { ci = d[ki][it[i-1]]+d[ki][it[i+1]] - d[it[i-1]][it[i+1]]; }
        /* Compute best edge "(it[jMin-1],it[jMin])" */
        /* and the cost "cjMin" of inserting "it[i]" into it: */
        { double da, db, cj; 
          long ja, jb, ka, kb;
          /* Consider inserting "it[i]" before "it[0]": */
          ja = (i == 0 ? 1 : 0);
          ka = it[ja]; 
          da = d[ki][ka];
          cjMin = da; jMin = ja;
          for (jb=ja+1; jb<N; jb++)
            { /* skip over "it[i]": */
              if (jb == i) continue;
              /* At this point "ka" is "it[ja]" where "ja" is */
              /* the previous value of "jb", skipping "i". */
              /* Consider inserting "it[i]" between "it[ja]" and "it[jb]": */
              kb = it[jb]; db = d[ki][kb];
              cj = da + db - d[ka][kb];
              if (cj < cjMin) { cjMin = cj; jMin = jb; }
              /* Prepare for next iteration: */
              ka = kb; da = db; ja = jb;
            }
          /* Consider inserting "it[i]" after end of path: */
          if (da < cjMin) { cjMin = da; jMin = N; }
        }
        /* Now dispose of "it[i]": */
        if (jMin == i) 
          { Error("program error - jMin=i"); }
        else if (ci > cjMin)
          { /* Move "ki = it[i]" to between "it[jMin-1" and "it[jMin]": */
            /* Don't increment "i", because "it[i]" changed. */
            if (jMin < i)
              { long r = i;
                while (r > jMin) { it[r] = it[r-1]; r--; }
                it[jMin] = ki;
              }
            else if (jMin > i+1)
              { long r = i+1;
                while (r < jMin) { it[r-1] = it[r]; r++; }
                it[jMin-1] = ki;
              }
            else
              { Error("program error - jMin = i+1, cjMin < ci"); }
          }
        else
          { i++; }
      }
  }
  
