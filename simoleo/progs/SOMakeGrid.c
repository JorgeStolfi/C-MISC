/* SOMakeGrid -- Create dyadic grids for testing. */
/* Last edited on 2006-03-19 12:44:58 by stolfi */

#include <SOGrid.h>
#include <SOParams.h>
#include <SOBasic.h>

#include <dg_grid.h>

#include <affirm.h>
#include <nat.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Creates grids for testing purposes.  The user may select
  between a few fixed geometries (a regular grid, a two-region
  one, etc.) or a random tree of a given maximum depth. */

/* INTERNAL PROTOTYPES */

typedef struct Options     /* Parsed command line options */
  { char *outName;
    dg_dim_t dim;        /* Grid dimension. */
    /* These are mutually exclusive: */
    bool_t regular;          /* TRUE asks for a regular grid. */
    bool_t circle;           /* TRUE for a grid that is denser within a circle. */
    bool_t front;            /* TRUE for a grid that is denser near a curved line. */
    int maxRank;           /* Maximum rank of grid. */
  } Options;

Options *GetOptions(int argn, char **argc);

SOGrid_Tree *SOMakeGrid_Make
  ( dg_dim_t dim,        /* Grid dimension. */
    bool_t regular,          /* TRUE asks for a regular grid. */
    bool_t circle,           /* TRUE for a grid that is denser within a circle. */
    bool_t front,            /* TRUE for a grid that is denser near a curved line. */
    int maxRank            /* Maximum rank of grid. */
  );
  
typedef bool_t SplitCriterion(dg_dim_t d, dg_rank_t r, double *lo, double *hi);
/* A predicate that tells whether the cell with rank {r},
  low corner {lo[0..d-1]}, and high corner {hi[0..d-1]} 
  should be subdivided. */

void SOMakeGrid_shatter_node
  ( dg_dim_t d,
    dg_tree_node_t *p,
    dg_rank_t r,
    double *lo, 
    double *hi,
    SplitCriterion split
  );

/* IMPLEMENTATIONS */

int main(int argn, char **argc)
  {
    Options *o = GetOptions(argn, argc);
    SOGrid_Tree *tree = 
      SOMakeGrid_Make(o->dim, o->regular, o->circle, o->front, o->maxRank);
    FILE *wr = open_write(txtcat(o->outName, ".tree"), TRUE);
    SOGrid_Tree_write(wr, tree);
    fclose(wr);
    return 0;
  }
    
SOGrid_Tree *SOMakeGrid_Make
  ( dg_dim_t dim,        /* Grid dimension. */
    bool_t regular,          /* TRUE asks for a regular grid. */
    bool_t circle,           /* TRUE for a grid that is denser within a circle. */
    bool_t front,            /* TRUE for a grid that is denser near a curved line. */
    int maxRank            /* Maximum rank of grid. */
  )
  { SOGrid_Tree *tree = SOGrid_Tree_new(dim);
    SplitCriterion *toobig;
    double minSize = pow(0.5, ((double)maxRank)/((double)dim)); // Min cell size.
    double lo[dim], hi[dim];
    int i;
    
    auto SplitCriterion regular_criterion, circle_criterion, front_criterion;

    bool_t regular_criterion(dg_dim_t d, dg_rank_t r, double *lo, double *hi)
      { return (r < maxRank); }
  
    bool_t circle_criterion(dg_dim_t d, dg_rank_t r, double *lo, double *hi)
      { double d2 = 0;  // Distance from root cell center to cell center, squared
        double sz = 0;  // Maximum extent of cell.
        int i;
	if (r >= maxRank) { return FALSE; }
        // Compute {d2}, {sz}:
        for (i = 0; i < d; i++)
          { double ci = (lo[i] + hi[i])/2.0 - 0.5;
            double szi = hi[i] - lo[i];
            d2 += ci*ci;
            if (szi > sz) { sz = szi; }
          }
        // Decide if cell is small enough:
        return (sz > (double)d2 * 0.8);
      }
  
    bool_t front_criterion(dg_dim_t d, dg_rank_t r, double *lo, double *hi)
      { double d2 = 0;  // Distance from root cell center to cell center, squared
        double sz = 0;  // Maximum extent of cell. 
        int i;
        if (r >= maxRank) { return FALSE; }
        // Compute {d2}, {sz}:
        for (i = 0; i < d; i++)
          { double ci = (lo[i] + hi[i])/2.0 - 0.5;
            double szi = hi[i] - lo[i];
            d2 += ci*ci;
            if (szi > sz) { sz = szi; }
          }
        // Decide if cell is small enough:
        { double dc2 = (d2 - 0.25)*(d2 - 0.25);
          return (sz > minSize + dc2/0.0625*(1.0 - minSize));
        }
      }
  
    if (regular)
      { toobig = &regular_criterion; }
    else if (circle)
      { toobig = &circle_criterion; }
    else if (front)
      { toobig = &front_criterion; }
    else
      { affirm(FALSE, "must specify some split criterion"); }
      
    for (i = 0; i < dim; i++) { lo[i] = 0.0; hi[i] = 1.0; }

    SOMakeGrid_shatter_node(dim, tree->root, 0, lo, hi, toobig);
    
    return tree;
  }

void SOMakeGrid_shatter_node
  ( dg_dim_t d,
    dg_tree_node_t *p,
    dg_rank_t r,
    double *lo, 
    double *hi,
    SplitCriterion *split
  )
  /*
    Given a node with no children, replaces it
    by a random subtree of depth at most maxRank. */
  {
    if (split(d, r, lo, hi))
      { double lo_child[d], hi_child[d];
        int i;
        dg_tree_node_split(p);
        // Compute low and high coordinates of low child:
        for (i = 0; i < d; i++)
          { lo_child[i] = lo[i]; 
            hi_child[i] = (i == (r%d) ? (lo[i]+hi[i])/2 : hi[i]);
          }
        SOMakeGrid_shatter_node(d, p->c[LO], r+1, lo_child, hi_child, split);
        
        // Compute low and high coordinates of high child:
        for (i = 0; i < d; i++)
          { lo_child[i] = (i == (r%d) ? (lo[i]+hi[i])/2 : lo[i]); 
            hi_child[i] = hi[i];
          }
        SOMakeGrid_shatter_node(d, p->c[HI], r+1, lo_child, hi_child, split);
      }
  }

#define PPUSAGE SOParams_SetUsage

Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);

    PPUSAGE(pp, "SOMakeGrid \\\n");
    PPUSAGE(pp, "  -outName NAME \\\n");
    PPUSAGE(pp, "  -dim D \\\n");
    PPUSAGE(pp, "  -maxRank NUM \\\n");
    PPUSAGE(pp, "  { -regular | -circle | -front }");

    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_GetKeyword(pp, "-dim");                               
    o->dim = SOParams_GetNextInt(pp, 1, 4);  

    SOParams_GetKeyword(pp, "-maxRank");                               
    o->maxRank = SOParams_GetNextInt(pp, 0, SOGRID_MAX_RANK);  

    o->regular = FALSE;
    o->circle = FALSE;
    o->front = FALSE;
    if (SOParams_KeywordPresent(pp, "-regular"))
      { o->regular = TRUE; }
    else if (SOParams_KeywordPresent(pp, "-circle"))
      { o->circle = TRUE; }
    else if (SOParams_KeywordPresent(pp, "-front"))
      { o->front = TRUE; }
    else
      { SOParams_Error(pp, "must specify type of grid"); }
    
    SOParams_Finish(pp);
    return o;
  }
