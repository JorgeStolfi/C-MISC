/* SOFindTentBasis -- Find a tent basis for a dyadic grid. */
/* Last edited on 2007-01-04 00:22:02 by stolfi */

#include <SOGrid.h>
#include <SOFunction.h>
#include <SOTentFunction.h>
#include <SOParams.h>
#include <SOBasic.h>

#include <dg_grid.h>

#include <stdio.h>
#include <values.h>

/* This program reads a grid {T} and computes a tent basis for the
  space of C0 degree 1 piecewise polynomials over that grid. */
  
typedef struct Options
  { char *treeName;   /* Name of grid tree file (minus the {.tree} extension). */
    char *outName;    /* Name of basis output file (minus the ".bas" extension) */
    bool_t maximal;     /* TRUE for maximal-support tents, FALSE for minimal ones. */
  } Options;
  
/* INTERNAL PROTOTYPES */

Basis SOFindTentBasis_make_basis(dg_dim_t d, SOTent_vec_t tents);

Options *GetOptions(int argn, char **argc);  

/* IMPLEMENTATIONS */

int main(int argn, char **argc)
  {
    Options *o = GetOptions(argn, argc);
    char *treeFile = txtcat(o->treeName, ".tree");
    FILE *rd = open_read(treeFile, TRUE);
    SOGrid_Tree *tree = SOGrid_Tree_read(rd);
    SOTent_vec_t tents;

    if (o->maximal)
      { tents = SOTent_maximal_basis(tree); }
    else
      { tents = SOTent_minimal_basis(tree); }

    { Basis bas = SOFindTentBasis_make_basis(tree->d, tents);
      FILE *wr = open_write(txtcat(o->outName, ".bas"), TRUE);
      SOFunction_WriteBasis(wr, bas);
      fclose(wr);
    }
    return 0;
  }
  
Basis SOFindTentBasis_make_basis(dg_dim_t d, SOTent_vec_t tents)
  { Basis bas = SOFunctionRef_vec_new(tents.nel);
    int i;
    for (i = 0; i < tents.nel; i++)
      { SOTent t = tents.el[i];
        bas.el[i] = (SOFunction *)SOTentFunction_FromHighBrick(d, t.cx);
      }
    return bas;
  }

#define PPUSAGE SOParams_SetUsage

Options *GetOptions(int argn, char **argc)
  {
    Options *o = (Options *)notnull(malloc(sizeof(Options)), "no mem");
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);
    PPUSAGE(pp, "SOFindTentBasis \\\n");
    PPUSAGE(pp, "  -treeName TRINAME \\\n");
    PPUSAGE(pp, "  [ -maximal | -minimal ] \\\n");
    PPUSAGE(pp, "  -outName BASISNAME\n");

    SOParams_GetKeyword(pp, "-treeName");                               
    o->treeName = SOParams_GetNext(pp);  

    if (SOParams_KeywordPresent(pp, "-maximal"))
      { o->maximal = TRUE; }
    else if (SOParams_KeywordPresent(pp, "-minimal"))
      { o->maximal = FALSE; }
    else
      { SOParams_Error(pp, "must specify \"-maximal\" or \"-minimal\""); }

    SOParams_GetKeyword(pp, "-outName");                               
    o->outName = SOParams_GetNext(pp);  

    SOParams_Finish(pp);
    return o;
  }
