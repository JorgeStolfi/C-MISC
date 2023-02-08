/* See SOIntegral.h */
/* Last edited on 2004-06-19 18:22:54 by stolfi */

#include <SOIntegral.h>
#include <SOBasic.h>

#include <dg_grid.h>

#include <affirm.h>
#include <nat.h>
#include <math.h>
#include <stdlib.h>

int SOIntegral_GaussOrder = 6;
  /* Number of knots to be used by {SOIntegral_Gauss} below,
    along each coordinate axis. */

/* Knots and weights for 2-point Gaussian quadrature: */
static double qpt_2[2] = 
  { -0.577350269189626, 
     0.577350269189626
  };
static double qwt_2[1] = 
  {  1.000000000000000 
  }; 

/* Knots and weights for 4-point Gaussian quadrature: */
static double qpt_4[4] = 
  { -0.861136311594053, 
     0.861136311594053, 
    -0.339981043584856, 
     0.339981043584856
  };
static double qwt_4[2] =
  { 0.347854845137454, 
    0.652145154862546
  };

/* Knots and weights for 6-point Gaussian quadrature: */
static double qpt_6[6] =  
  { -0.932469514203152, 
     0.932469514203152, 
    -0.661209386466265,
     0.661209386466265, 
    -0.238619186083197, 
     0.238619186083197
  };
static double qwt_6[3] = 
  { 0.171324492379170, 
    0.360761573048139, 
    0.467913934572691 
  };

/* Knots and weights for 8-point Gaussian quadrature: */
static double qpt_8[8] =
  { -0.960289856497536,
     0.960289856497536,
    -0.796666477413627,
     0.796666477413627,
    -0.525532409916329,
     0.525532409916329,
    -0.183434642495650,
     0.183434642495650
  };
static double qwt_8[4] = 
  { 0.101228536290376,
    0.222381034453374,
    0.313706645877887,
    0.362683783378362
  };

/* Knots and weights for 10-point Gaussian quadrature: */
static double qpt_10[10] =
  { -0.973906528517172,
     0.973906528517172,
    -0.865063366688985,    
     0.865063366688985, 
    -0.679409568299024,
     0.679409568299024,
    -0.433395394129247,
     0.433395394129247,
    -0.148874338981631,
     0.148874338981631
  };
static double qwt_10[5] = 
  { 0.066671344308688,
    0.149451349150581,
    0.219086362515982, 
    0.269266719309996,
    0.295524224714753
  };

/* Knots and weights for 12-point Gaussian quadrature: */
/*static double qpt_12[12] =
  { -0.981 560 634 246 719,
     0.981 560 634 246 719,
    -0.904 117 256 370 475,
     0.904 117 256 370 475,
    -0.769 902 674 194 305,
     0.769 902 674 194 305,
    -0.587 317 954 286 617,
     0.587 317 954 286 617,
    -0.367 831 498 998 180,
     0.367 831 498 998 180,
    -0.125 233 408 511 469,
     0.125 233 408 511 469
  };
static double qwt_12[6] = 
  { 0.047 175 336 386 512,
    0.106 939 325 995 318,
    0.160 078 328 543 346,
    0.203 167 426 723 066,
    0.233 492 536 538 355,
    0.249 147 045 813 403
  };
*/

/* Pointers to knot and weight tables: */
static double *qpoints[5] = 
  { (double *)qpt_2,
    (double *)qpt_4,
    (double *)qpt_6,
    (double *)qpt_8,
    (double *)qpt_10
  };
static double *qweight[5] =
  { (double *)qwt_2,
    (double *)qwt_4,
    (double *)qwt_6,
    (double *)qwt_8,
    (double *)qwt_10
  };


void SOIntegral_Gauss
  ( SOIntegral_Func f, 
    dg_dim_t d, 
    dg_dim_t fd, 
    double *sum, 
    double *corr
  )
  { int ktable = (SOIntegral_GaussOrder/2)-1;   /* Which tables to use */
    double *pt = qpoints[ktable]; /* Knot positions in [-1_+1]. */
    double *wt = qweight[ktable];
    double volfactor = (double)1.0/((double)(1<<d)); /* {2^{-d}}: */
    static int ix[MAX_PDIM];    /* Counters for the grid. */
    static double x[MAX_PDIM];  /* Coordinates of sample points, in [0_1]. */
    static double fx[MAX_FDIM]; /* Function value at {x}. */
    int i, j;  
    
    /* Initialize counter variables: */
    for (i = 0; i < d; i++) { ix[i] = 0;  x[i] = (double)(1 + pt[0])/2; }
    
    /* Simulate {d} nested loops, with {ix[i]} varying from 0 to {SOIntegral_GaussOrder-1}: */
    do 
      { 
        /* Compute weight {w} of sample point, including the box measure: */
        double w = volfactor;
	for (i = 0; i < d; i++) { w = w * wt[ix[i]/2]; }
       
        /* Evaluate the function at the sample point: */
        f(&(x[0]), &(fx[0]));
        
        /* Accumulate {w*f(x)} on the integral: */
         for (j = 0; j < fd; j++)
          { double term = w * fx[j];
            /* Kahan's summation formula: */
            double sumj = sum[j];
            double tcorr = term - corr[j];
            double newSum = sumj + tcorr;
            corr[j] = (newSum - sumj) - tcorr;
            sum[j] = newSum;
          }
        
        /* Increment the loop counters and update sample point: */
        i = d;
        do
          { i--;
            ix[i] = (ix[i] + 1) % SOIntegral_GaussOrder;
            x[i] = (double)(1 + pt[ix[i]])/2; 
          }
        while ((i > 0) && (ix[i] == 0));
        /* After the last point, {i = ix[0] = 0}. */
      }
    while ((i > 0) || (ix[i] != 0));  

  }


void SOIntegral_Gauss_Limits
  ( SOIntegral_Func f, 
    dg_dim_t d, 
    dg_dim_t fd, 
    double *sum, 
    double *corr,
    double *min,    /* Lower coordinates of integration limits. */
    double *max     /* Higher coordinates of integration limits. */

  )
  { int ktable = (SOIntegral_GaussOrder/2)-1;   /* Which tables to use */
    double *pt = qpoints[ktable]; /* Knot positions in [-1_+1]. */
    double *wt = qweight[ktable];
    double volfactor = 1;
    static int ix[MAX_PDIM];    /* Counters for the grid. */
    static double x[MAX_PDIM];  /* Coordinates of sample points, in [0_1]. */
    static double fx[MAX_FDIM]; /* Function value at {x}. */
    int i, j;  
    
    /* Initialize counter variables: */
    for (i = 0; i < d; i++) 
      { 
        volfactor *= (double)(max[i]-min[i]) / 2;
        ix[i] = 0;  
        x[i] = (double)((max[i]+min[i])+(double)(max[i]-min[i])*pt[0])/2; 
      }
    
    /* Simulate {d} nested loops, with {ix[i]} varying from 0 to {SOIntegral_GaussOrder-1}: */
    do 
      { 
        /* Compute weight {w} of sample point, including the box measure: */
        double w = volfactor;
	for (i = 0; i < d; i++) { w = (double)w * wt[ix[i]/2]; }
       
        /* Evaluate the function at the sample point: */
        f(&(x[0]), &(fx[0]));
        
        /* Accumulate {w*f(x)} on the integral: */
         for (j = 0; j < fd; j++)
          { double term = (double)w * fx[j];
            /* Kahan's summation formula: */
            double sumj = sum[j];
            double tcorr = term - corr[j];
            double newSum = sumj + tcorr;
            corr[j] = (newSum - sumj) - tcorr;
            sum[j] = newSum;
          }
        
        /* Increment the loop counters and update sample point: */
        i = d;
        do
          { i--;
            ix[i] = (ix[i] + 1) % SOIntegral_GaussOrder;
            x[i] = (double)((max[i]+min[i])+(double)(max[i]-min[i])*pt[ix[i]])/2;
          }
        while ((i > 0) && (ix[i] == 0));
        /* After the last point, {i = ix[0] = 0}. */
      }
    while ((i > 0) || (ix[i] != 0));  

  }


void SOIntegral_OnRootCell
  ( SOIntegral_Func f, 
    dg_dim_t d,       /* Dimension {m} of domain. */
    dg_dim_t fd,      /* Dimension {n} of the function's range. */
    SOGrid_Tree *tree,  /* Finite cell grid to use. */
    double *sum, 
    double *corr
  )
  {
    /* If single cell, use SOIntegral_Gauss, else recurse on
       each sub-cell and add results. */
    //    if (tree == NULL)
    //      { SOIntegral_Gauss(f, d, fd, sum, corr); }
    //    else
    { affirm(FALSE, "not implemented yet"); }
  }

