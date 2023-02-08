/* See {output.h}. */
/* Last edited on 2009-02-10 10:05:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>

#include <output.h>

#include <r3.h>
#include <jsfile.h>
#include <jsstring.h>

#include <basic.h>
#include <geomodel.h>
#include <wavefront.h>
#include <mesh.h>

/* INTERNAL PROTOTYPES */
char *iteration_tag(int iter);

void plot_wave(wavefront_t *wf, char *outName, int iter)
  { }  
  
/* IMPLEMENTATIONS */

void output_wave(wavefront_t *wf, char *outName, int iter)
  {
    /* Open output file: */
    char *iterTag = iteration_tag(iter);
    char *fname = NULL;
    asprintf(&fname, "%s-%s.tri", outName, iterTag);
    FILE *wr = open_write(fname, TRUE);
    
    /* Set sample numbers: */ 
    int i;
    for (i = 0; i < wf->ns; i++)
      { segment_t *si = wf->st.e[i];
        si->num = i;
      }
    
    /* Write wavefront mesh to file {wr}: */ 
    mesh_t *tri = mesh_from_topology(wf->a);
    write_mesh(wr, tri);
    
    /* Close file: */
    fclose(wr);
    free(fname);
    free(iterTag);
  } 

void output_model(geomodel_t *geo, char *outName)
  {
    /* Open output file: */
    char *fname = NULL;
    asprintf(&fname, "%s.geo", outName);
    FILE *wr = open_write(fname, TRUE);
    
    /* Write geophysical model to file {wr}: */ 
    write_geomodel(wr, geo);
    
    /* Close file: */
    fclose(wr);
    free(fname);
  } 
  
char *iteration_tag(int iter)
  { if (iter == -1) 
      { return txtcat("ini",""); }
    else if (iter == INT_MAX)
      { return txtcat("fin",""); }
    else 
      { affirm(iter > 0, "bad iteration number"); 
        return fmtint(iter,6);
      }
  }

