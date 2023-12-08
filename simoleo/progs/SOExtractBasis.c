/* SOExtractBasis -- Extracts elements from a basis */
/* Last edited on 2023-10-15 03:38:48 by stolfi */

#include <SOParams.h>
#include <SOFunction.h> 
#include <SOSpline.h> 
#include <SOBasic.h>
#include <vec.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <stdio.h>

typedef struct Options
  { char *baseName;     /* Name of basis file (minus ".bas" extension). */
    char *outName;      /* Prefix for output file names. */
    nat_vec_t index;    /* Indices of elements to extract. */
    nat_vec_t vertex;   /* Extract elements associated with these vertices. */
    nat_vec_t edge;     /* Extract elements associated with these edges. */
    nat_vec_t face;     /* Extract elements associated with these faces. */
    bool_t writeElems;    /* TRUE to write extracted elements as individual files. */
    bool_t writeBasis;    /* TRUE to write extracted elements as a basis. */
    bool_t writeRest;     /* TRUE to write basis minus extracted elements. */
  } Options;

Options GetOptions(int argn, char **argc);

void CompressBasis(Basis *B);
  /* Compress non-NULL elements of {B} to initial positions, and
    trims it. */

Basis ReadBasis(char *name);
  /* Reads a function basis from the file {name} plus extension ".bas". */
  
void WriteBasis(Basis F, char *name);
  /* Writes basis {F} to the file {name} plus extension ".bas". */
  
void WriteFunction(SOFunction *f, char *name);
  /* Writes function {f} to the file {name} plus extension ".fun". */
  
bool_t Selected(SOFunction *f, nat_vec_t vertex, nat_vec_t edge, nat_vec_t face);
  /* TRUE if {f} is an instance of {SOSpline}, and 
    the topological elements (vertices, edges,or faces) shared by all 
    pieces of {f} include any of the elements listed in {vertex}, 
    {edge}, or {face}. */
  
int main(int argn, char **argc)
  { Options o = GetOptions(argn, argc);
    FILE *rd = open_read(txtcat(o.baseName, ".bas"), TRUE);
    Basis F = SOFunction_ReadBasis(rd);
    Basis G = SOFunctionRef_vec_new(F.ne);
    int i, k;
    for (k = 0; k < G.ne; k++) { G.e[k] = NULL; }
    /* Extract elements specified by their indices: */
    for (i = 0; i < o.index.ne; i++)
      { int k = o.index.e[i];
        G.e[k] = F.e[k];
        F.e[k] = NULL; /* For now... */
      }
    /* Extract elements specified by site/edge/face number: */
    for (k = 0; k < F.ne; k++)
      { if ((F.e[k] != NULL) && Selected(F.e[k], o.vertex, o.edge, o.face))
          { G.e[k] = F.e[k];
            F.e[k] = NULL;  /* For now... */
          }
      }
    if (o.writeElems)
      { for (k = 0; k < G.ne; k++)
          { if (G.e[k] != NULL)
              { char *elName = txtcat3(o.outName, "-e", fmt_int(k, 4));
                WriteFunction(G.e[k], elName);
              }
          }
      }
    if (o.writeRest)
      { CompressBasis(&F);
        WriteBasis(F, txtcat(o.outName, "-r"));
      }
    if (o.writeBasis)
      { CompressBasis(&G);
        WriteBasis(G, txtcat(o.outName, "-x"));
      }
    return 0;
  }
  
void CompressBasis(Basis *B)
  { int k, newN = 0;
    for (k = 0; k < B->ne; k++)
      { if (B->e[k] != NULL) { B->e[newN] = B->e[k]; newN++; } }
    SOFunctionRef_vec_trim(B, newN);
  }
  
Basis ReadBasis(char *name)
  { FILE *rd = open_read(txtcat(name, ".bas"), TRUE);
    Basis F = SOFunction_ReadBasis(rd);
    fclose(rd);
    return F;
  }
  
void WriteBasis(Basis F, char *name)
  { FILE *wr = open_write(txtcat(name, ".bas"), TRUE);
    SOFunction_WriteBasis(wr, F);
    fclose(wr);
  }
  
void WriteFunction(SOFunction *f, char *name)
  { FILE *wr = open_write(txtcat(name, ".fun"), TRUE);
    f->m->write(f, wr);
    fclose(wr);
  }

bool_t Selected(SOFunction *f, nat_vec_t vertex, nat_vec_t edge, nat_vec_t face)
  { SOSpline *fpw = SOSpline_Cast(f);
    if (fpw != NULL)
      { SOGrid *tree = fpw->d->tree;
        PieceDataRef_vec_t *pd = &(fpw->d->pd);
        SiteNumber vCom[3]; int nvCom;
        EdgeNumber eCom[3]; int neCom;
        FaceNumber fCom[3]; int nfCom;
        int i, j;
        FindSharedElems(pd, tree, vCom, &nvCom, eCom, &neCom, fCom, &nfCom);
        if (nvCom == 1)
          { for (i = 0; i < nvCom; i++)
              for (j = 0; j < vertex.ne; j++)
                { if (vCom[i] == vertex.e[j]) { return TRUE; } }
          }
        if (neCom == 1)
          { for (i = 0; i < neCom; i++)
              for (j = 0; j < edge.ne; j++)
                { if (eCom[i] == edge.e[j]) { return TRUE; } }
          }
        if (nfCom == 1)
          { for (i = 0; i < nfCom; i++)
              for (j = 0; j < face.ne; j++)
                { if (fCom[i] == face.e[j]) { return TRUE; } }
          }
      }
    return FALSE;
  }    

#define PPUSAGE SOParams_SetUsage

nat_vec_t ParseNumList(SOParams_T *pp, char *key);
  /* Parses the command line arguments for zero or more instances 
    of {key} followed by an integer, and returns a vector of those
    integers. */

Options GetOptions(int argn, char **argc)
  {
    Options o;
    SOParams_T *pp = SOParams_NewT(stderr, argn, argc);
          
    PPUSAGE(pp, "SOExtractBasis \\\n");
    PPUSAGE(pp, "  -baseName NAME \\\n");
    PPUSAGE(pp, "  -outName NAME \\\n");
    PPUSAGE(pp, "  [-writeElems] [-writeBasis] [-writeRest] \\\n");
    PPUSAGE(pp, "  [ -index NUM ]...\n");

    SOParams_GetKeyword(pp, "-baseName");                               
    o.baseName = SOParams_GetNext(pp);  

    SOParams_GetKeyword(pp, "-outName");                               
    o.outName = SOParams_GetNext(pp);  

    o.writeElems = SOParams_KeywordPresent(pp, "-writeElems");
    o.writeBasis = SOParams_KeywordPresent(pp, "-writeBasis");
    o.writeRest = SOParams_KeywordPresent(pp, "-writeRest");
    if (! (o.writeElems || o.writeBasis || o.writeRest))
      { SOParams_Error(pp, 
          "must specify \"-writeElems\", \"-writeBasis\", or \"-writeRest\""
        );
      }
    
    o.index = ParseNumList(pp, "-index");
    o.vertex = ParseNumList(pp, "-vertex");
    o.edge = ParseNumList(pp, "-edge");
    o.face = ParseNumList(pp, "-face");
    
    SOParams_Finish(pp);
    return o;
  }

nat_vec_t ParseNumList(SOParams_T *pp, char *key)
  { nat_vec_t num = nat_vec_new(10);
    int nNums = 0;
    while (SOParams_KeywordPresent(pp, key))
      { nat_vec_expand(&num, nNums);
        num.e[nNums] = SOParams_GetNextInt(pp, 0, INT_MAX);
        nNums++;
      }
    nat_vec_trim(&num, nNums);
    return num;
  }

                        
