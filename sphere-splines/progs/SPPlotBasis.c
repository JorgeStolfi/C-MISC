/* SPPlotBasis.c -- Plots elements from a spherical function basis. */
/* Last edited on 2008-05-24 12:28:19 by stolfi */

#include <SPBasic.h>
#include <vec.h>
#include <SPOptions.h>
#include <SPFunction.h>
#include <SPSpline.h>
#include <SPTriang.h>
#include <SPIntegral.h>
#include <SPPlot.h>
#include <SPPlotOptions.h>
#include <SPH3.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <pswr.h>
#include <r3.h>
#include <r4.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <values.h>
#include <limits.h>

typedef struct Options /* Parsed command line options */
  { char *basisName;      /* Name of basis file (minus ".bas" extension). */
    bool_t showTriang;      /* If elem is an {SPSpline}, show its triang.  */
    bool_t showSupp;        /* If elem is an {SPSpline}, show supp. plane.  */
    char *outName;        /* Prefix for output file names. */
    /* Element selection parameters: */
    bool_t all;             /* TRUE to plot all elements. */
    int_vec_t index;        /* Indices of elements to extract. */
    int_vec_t face;         /* Extract elements associated with these faces. */
    int_vec_t edge;         /* Extract elements associated with these edges. */
    int_vec_t vertex;       /* Extract elements associated with these vertices. */
    /* Plotting options: */
    SPPlotOptions_t plt;    /* Plot format and style options. */
  } Options;

Options GetOptions(int argn, char **argc);

Basis ReadBasis(char *name);
  /* Reads a function basis from the file {name} plus extension ".bas". */
  
  /* TRUE if {f} is an instance of {SPSpline}, and 
    the topological elements (vertices, edges,or faces) shared by all 
    pieces of {f} include any of the elements listed in {vertex}, 
    {edge}, or {face}. */
 
Triangulation *GetTriangulation(SPFunction *f);
SPH3_Plane GetSupportingPlane(SPFunction *f);

void ExtractAllElems
  ( Basis F, 
    int_vec_t *idX, int *nX
  );
    
void ExtractSpecificElems
  ( Basis F, 
    int_vec_t index, 
    int_vec_t *idX, int *nX
  );
  
void ExtractElemsBySupport
  ( Basis F, 
    Triangulation *tri,
    int dim, 
    int_vec_t num, 
    int_vec_t *idX, int *nX
  );
  
bool_t SupportSelected
  ( PieceDataRef_vec_t *pd,
    Triangulation *tri,
    int dim,
    int num
  );
  /* True if a basis element with pieces {pd} is associated to
    a topological element of dimension {dim} whose number is {num}.
    Specifically, if

      {dim = 2} and {pd} is the single face {num}; or
      {dim = 1} and {pd} is two faces sharing edge {num}; or
      {dim = 0} and {pd} is 3 or more faces sharing vertex {num}.

  */
    
void SelectCommonScale
  ( Basis F, 
    int_vec_t idX, 
    double *fRange, 
    double *fStep
  ); 
  /* Selects common values of {fRange} and {fStep} that are
     adequate for plotting all of the members of {F} listed
     in {idX}. */

void PlotElement
  ( SPPlot_Stream *fps,
    int index,
    SPFunction *f,
    Options *o 
  );
  /* Plots element {f} of basis, whose index is {index}. */

int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    SPPlotOptions_t *po = &(o.plt);
    int smpOrder = 8; /* Rough integrals are OK here. */
    SPIntegral_SetDefaultSamplingOrder(smpOrder);
    FILE *rd = open_read(txtcat(o.basisName, ".bas"), TRUE);
    Basis F = SPFunction_ReadBasis(rd);
    Triangulation *tri = SPSpline_BasisTriangulation(F);
    int_vec_t idX = int_vec_new(10);
    int i;

    /* Put indices of selected elements in {idX} array, in plotting order. */
    /* Leave an invalid index between each group, as page break marks. */
    { int nX = 0;
      if (o.all) 
        { ExtractAllElems(F, &idX, &nX); }
      else
        { ExtractSpecificElems(F, o.index, &idX, &nX);
          ExtractElemsBySupport(F, tri, 2, o.face,   &idX, &nX);
          ExtractElemsBySupport(F, tri, 1, o.edge,   &idX, &nX);
          ExtractElemsBySupport(F, tri, 0, o.vertex, &idX, &nX);
        }
      int_vec_trim(&idX, nX);
    }
    
    /* Fix global isoline parameters, if needed: */
    if (po->fRange <= 0.0)
      { SelectCommonScale(F, idX, &(po->fRange), &(po->fStep)); }

    /* Fix global perspective parameters, if needed: */
    SPPlot_FixView(&(po->obs), 0.0, &(po->upp), &(po->radius), 0.0, NULL, TRUE);
    
    /* Fix figure size: */
    if (po->figSize == 0.0)
      { po->figSize = (po->eps ? 150.0 : 70.0); }
      
    /* Open the figure stream: */
    SPPlot_Stream *fps = SPPlot_NewStream
      (po->eps, o.outName, po->paperSize, po->figSize, po->projLonLat, po->caption.ne);

    /* Now plot the selected elements: */
    for (i = 0; i < idX.ne; i++)
      { int k = idX.e[i];
        if ((k >= 0) && (k < F.ne)) 
          { SPFunction *Fk = F.e[k];
            if (Fk != NULL)
              { PlotElement(fps, k, Fk, &o);
                if (po->eps) { /* Avoid duplicates: */ F.e[k] = NULL; }
              }
          }
        else
          { if (! po->eps) 
              { /* Finish the current row in page: */
                pswr_fill_row(fps);
              }
          }
      }
    pswr_close_stream(fps);
    return 0;
  }
  
Basis ReadBasis(char *name)
  { FILE *rd = open_read(txtcat(name, ".bas"), TRUE);
    Basis F = SPFunction_ReadBasis(rd);
    fclose(rd);
    return F;
  }
  
void ExtractAllElems
  ( Basis F, 
    int_vec_t *idX, int *nX
  )
  { int k;
    for (k = 0; k < F.ne; k++)
      { if (F.e[k] != NULL) 
          { int_vec_expand(idX, *nX);
            idX->e[*nX] = k; (*nX)++;
          }
      }
    int_vec_expand(idX, *nX);
    idX->e[*nX] = F.ne; (*nX)++;
  }

void ExtractSpecificElems
  ( Basis F, 
    int_vec_t index, 
    int_vec_t *idX, int *nX
  )
  { int i;
    for (i = 0; i < index.ne; i++)
      { int k = index.e[i];
        affirm((k >= 0) && (k < F.ne), "invalid element index"); 
        if (F.e[k] != NULL) 
          { int_vec_expand(idX, *nX);
            idX->e[*nX] = k; (*nX)++;
          }
      }
    int_vec_expand(idX, *nX);
    idX->e[*nX] = F.ne; (*nX)++;
  }

void ExtractElemsBySupport
  ( Basis F, 
    Triangulation *tri,
    int dim, 
    int_vec_t num, 
    int_vec_t *idX, int *nX
  )
  { int j, k;
    for (j = 0; j < num.ne; j++) 
      { int numj = num.e[j];
        for (k = 0; k < F.ne; k++)
          { SPSpline *Fk = SPSpline_Cast(F.e[k]);
            if 
              ( (Fk != NULL) && 
                (Fk->d->tri == tri) &&
                SupportSelected(&(Fk->d->pd), tri, dim, numj)
              )
              { 
                fprintf(stderr, "  bas[%04d] tri elem %04d dim = %d\n", k, numj, dim);
                int_vec_expand(idX, *nX);
                idX->e[*nX] = k; (*nX)++;
              }
          }
	fprintf(stderr, "  ---\n");
        int_vec_expand(idX, *nX);
        idX->e[*nX] = F.ne; (*nX)++;
      }
  }

bool_t SupportSelected
  ( PieceDataRef_vec_t *pd,
    Triangulation *tri,
    int dim,
    int num
  )
  { SiteNumber vCom[3]; int nvCom;
    EdgeNumber eCom[3]; int neCom;
    FaceNumber fCom[3]; int nfCom;
    int i;
    FindSharedElems(pd, tri, vCom, &nvCom, eCom, &neCom, fCom, &nfCom);
    switch (dim)
      { 
        case 0:
          if (nvCom != 1) { /* Not a vertex element */ return FALSE; }
          affirm((neCom == 0) && (nfCom == 0), "funny support");
          for (i = 0; i < nvCom; i++)
            if (vCom[i] == num) { return TRUE; }
          return FALSE;

        case 1:
          if (neCom != 1) { /* Not an edge element */ return FALSE; }
          affirm((nvCom == 2) && (nfCom == 0), "funny support");
          for (i = 0; i < neCom; i++)
            if (eCom[i] == num) { return TRUE; }
          return FALSE;
          
        case 2:
          if (nfCom != 1) { /* Not a face element */ return FALSE; }
          for (i = 0; i < nfCom; i++)
            if (fCom[i] == num) { return TRUE; }
          return FALSE;
          
        default:
          affirm(FALSE, "bad dim");
          return FALSE;
      }
  }   

void SelectCommonScale
  ( Basis F, 
    int_vec_t idX, 
    double *fRange, 
    double *fStep
  )
  { double fMax2Sum = 0.0;
    double fMaxMax = 0.0;
    int i;
    fprintf(stderr, "computing common range...\n");
    for (i = 0; i < idX.ne; i++)
      { int k = idX.e[i];
        if ((k >= 0) && (k < F.ne)) 
          { SPFunction *fk = F.e[k];
            if (fk != NULL) 
              { double fkMax = SPPlot_EstRange(fk, FALSE);
                if (fkMax > fMaxMax) { fMaxMax = fkMax; }
                fMax2Sum += fkMax*fkMax;
              }
          }
      }
    /* Pick maximum of all {fkMax}, but beware of way-out weirdos: */
    { double fMaxAvg = sqrt(fMax2Sum/F.ne); 
      double fLim = 2.0*fMaxAvg;
      double idealStep = (fMaxMax < fLim ? fMaxMax : fLim)/DefaultIsolines;
      (*fStep) = SPPlot_RoundToNice(idealStep);
      (*fRange) = DefaultIsolines * (*fStep);
    }
  }

void PlotElement
  ( SPPlot_Stream *fps,
    int index,
    SPFunction *f,
    Options *o
  )
  { Triangulation *tri = GetTriangulation(f);
    SPPlotOptions_t *po = &(o->plt);
    SPH3_Plane supp = (o->showSupp ? GetSupportingPlane(f) : Omega);
    SPH3_Point elObs = po->obs, elUpp = po->upp;
    double elRadius = po->radius;
    double fMinObs = INFINITY, fMaxObs = -INFINITY;
    double relMeshSize;

    auto double evalFunction(S2Point *p);

    double evalFunction(S2Point *p)
      { return f->m->eval(f, p); }
   
    /* Adjust observer and zenith, for this element: */
    SPPlot_FixView
      (&(elObs), po->autoObs, &(elUpp), &(elRadius), po->autoRadius, f, TRUE);
      
    /* Compute mesh size in world coordinates: */
    relMeshSize = elRadius * po->meshSize/(po->figSize/2);

    /* Fix function range, if needed: */
    SPPlot_FixRange(&(po->fRange), &(po->fStep), po->autoRange, f, TRUE);

    /* Set scales and start new figure: */
    double xMax = 1.0;
    double yMax = (po->projLonLat ? 0.5 : 1.0);
    pswr_new_picture(fps, -xMax,+xMax, -yMax,+yMax);
    
    /* Plot element: */
    SPPlot_Everything
      ( fps,
        evalFunction, 
        tri, 
        relMeshSize, 
        o->showTriang, 
        &supp,
        po->fRange, po->fStep,
        &(elObs), &(elUpp), 
        elRadius,
        &(po->light),
        po->lineWidth, 
        0, 0, FALSE,
        TRUE,
        &fMinObs, &fMaxObs
      );
    if (! po->eps)
      { 
        /* Make caption: */
        string_vec_t icap = SPPlot_ExpandCaption
          (&po->caption, +1, index, 0.0, 0.0, po->fRange, po->fStep, "");
        int k;
        for (k = 0; k < icap.ne; k++)
          { pswr_add_caption(fps, icap.e[k], 0.5); }
        SPPlot_FreeCaption(&icap);
      }
    fprintf(stderr, "observed range = [ %g _ %g ]\n", fMinObs, fMaxObs);
  }

Triangulation *GetTriangulation(SPFunction *f)
  { SPSpline *fpw;
    if ((fpw = SPSpline_Cast(f)) != NULL) 
      { return fpw->d->tri; }
    else
      { fprintf(stderr, "** warning - no triangulation - type = %s\n", f->type);
        return NULL;
      }
  }

SPH3_Plane GetSupportingPlane(SPFunction *f)
  { SPSpline *fpw;
    if ((fpw = SPSpline_Cast(f)) != NULL)
      { return fpw->d->supp; }
    else
      { fprintf(stderr, "** warning - no supporting plane - type = %s\n", f->type);
        return Omega;
      }
  }


Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);

    SPOptions_SetUsage(pp, 
      "SPPlotBasis \\\n"
      "  -basisName NAME \\\n"
      "  [ -showTriang ] [ -showSupp ] \\\n"
      "  -outName NAME \\\n"
      "  [ -all ] [ -index NUM ]... \\\n"
      "  [ -face NUM ]... [ -edge NUM ]... [ -vertex NUM ]... \\\n"
      SPPlotOptions_FunctionHelp " \n"
    );

    if (SPOptions_TestKeyword(pp, "-basisName"))
      { o.basisName = SPOptions_GetNext(pp); }
    else
      { o.basisName = "-"; }  

    o.showTriang = SPOptions_TestKeyword(pp, "-showTriang");
    
    o.showSupp = SPOptions_TestKeyword(pp, "-showSupp");
    
    SPOptions_GetKeyword(pp, "-outName");                               
    o.outName = SPOptions_GetNext(pp);
    
    /* Element selection parameters: */
    
    o.all = SPOptions_TestKeyword(pp, "-all");
    o.index = SPOptions_GetIntList(pp, "-index", 0, INT_MAX);
    o.face = SPOptions_GetIntList(pp, "-face", 0, INT_MAX);
    o.edge = SPOptions_GetIntList(pp, "-edge", 0, INT_MAX);
    o.vertex = SPOptions_GetIntList(pp, "-vertex", 0, INT_MAX);
    
    o.plt = SPPlotOptions_FunctionParse(pp);

    SPOptions_Finish(pp);
    return o;
  }
