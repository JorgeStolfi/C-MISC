/* Last edited on 2004-06-19 22:34:38 by stolfi */

----------------------------------------------------------------------
SOMakeWaveBasis.c 
  
FILE *OpenPlot(char *name)
  {
    double min[2], max[2];
    FILE *plot_file = open_write(name);
    min[X] = -1.05; max[X] = +1.05;
    min[Y] = -1.05; max[Y] = +1.05;
    SO2DPlot_Begin_Figure( plot_file, min, max );    
    eps_draw_segment(plot_file, min[X], 0.0, max[X], 0.0);
    eps_draw_segment(plot_file, 0.0, min[Y], 0.0, max[Y]);
    eps_draw_rectangle(plot_file, min[X], max[X], min[Y], max[Y]);
    return plot_file;
  }
    
----------------------------------------------------------------------

typedef struct SOGrid_Node
  { struct SOGrid_Node *p;    /* Parent node. */
    struct SOGrid_Node *c[2]; /* Children nodes (both NULL or both non-NULL). */
    dg_cell_index_t index;       /* Node index, as in the infinite binary tree. */
    dg_cell_index_t count;       /* Number of nodes in the subtree rooted at this. */
  } SOGrid_Node;
  /* A {SOGrid_Node} may represent, for instance, a cell of the dyadic cell
    tree. The children {c[LO],c[HI]} would then correspond to the two
    half-cells that result from bisecting the parent cell transversely
    to the axis of longest extent.  A childless node (a `leaf') represents 
    a cell that was left undivided. 
    
    In that interpretation, a non-leaf node is a super-cell which
    contains all its descendant cells, and is the union of all its
    descendant leaf cells.
    
    A {SOGrid_Node} may also be interpreted the point that is the center
    of that cell. (Note that there is at most one cell of the dyadic
    cell tree which is centered at a given a point {p} of {R^d}.) Then
    the leaf nodes, in this interpretation, correspond to the centers
    of undivided cells of the former interpretation; and internal
    nodes correspond to centers of divided cells. */

SOGrid_Node *SOGrid_Node_new(SOGrid_Node *p, dg_cell_index_t k);
  /* Creates a new node with no children, parent {p},
    and index {k}. */
  
void SOGrid_Node_free(SOGrid_Node *n);
  /* Reclaims the storage used by {n}, which must not have
    any children.  Also reclaims {n->data}, if present. */

void SOGrid_Node_split(SOGrid_Node *n);
  /* The node {n} must have no children. Creates two new nodes and
    appends them as children of {n}. Note that the {count} fields of
    {n}'s ancestors will need fixing after the split. */

void dg_tree_node_unsplit(SOGrid_Node *n);
  /* The node {n} must have two children but no grandchildren. Deletes
    and reclaims the two children of {n}, leaving {n} childless. Note
    that the {count} fields of {n}'s ancestors will need fixing after
    the split. */

void SOGrid_free_subtree(SOGrid_Node* n);
  /* Recursively reclaims {n} and its descendents. */

/* Vectors of {SOGrid_Node*}: */

typedef struct SOGrid_NodeRefVec { nat nel; SOGrid_Node **el; } SOGrid_NodeRefVec;

#define SOGrid_NodeRefVec_expand(nv,index) \
  vec_expand(vec_cast_ref(nv), index, sizeof(SOGrid_Node *))

#define SOGrid_NodeRefVec_trim(nv,nel) \
  vec_trim(vec_cast_ref(nv), nel, sizeof(SOGrid_Node *))


----------------------------------------------------------------------

void SOGrid_Node_split(SOGrid_Node *n)
  { affirm(n->c[0] == NULL, "node has c[LO]");
    affirm(n->c[1] == NULL, "node has c[HI]");
    n->c[0] = SOGrid_Node_new(n, 2 * n->index);
    n->c[1] = SOGrid_Node_new(n, 2 * n->index + 1);
    n->count = 3;
  }

void SOGrid_Node_unsplit(SOGrid_Node *n)
  { affirm(n->c[0] != NULL, "node has no c[LO]");
    affirm(n->c[1] != NULL, "node has no c[HI]");
    SOGrid_Node_free(n->c[0]); n->c[0] = NULL; 
    SOGrid_Node_free(n->c[1]); n->c[1] = NULL;
  }

  
void SOGrid_subtree_free(SOGrid_Node* root)
  { if (root != NULL)
      { SOGrid_subtree_free(root->c[0]);
        SOGrid_subtree_free(root->c[1]);
        SOGrid_Node_free(root);
      }
  }

SOGrid_Node *SOGrid_Node_new(SOGrid_Node *p, dg_cell_index_t k)
  { void *v = notnull(malloc(sizeof(SOGrid_Node)), "no mem for SOGrid_Node");
    SOGrid_Node *n = (SOGrid_Node*)v;
    n->p = p;
    n->c[0] = NULL;
    n->c[1] = NULL;
    n->index = k;
    n->count = 1;
    return n;
  }

void SOGrid_Node_free(SOGrid_Node *node)
  { affirm(node->c[0] == NULL, "node has c[LO]");
    affirm(node->c[1] == NULL, "node has c[HI]");
    free(node);
  }

/*
SOGrid_NodeRefVec SOGrid_NodeRefVec_new(nat nel)
  {  // This is not a macro only because gcc does not allow cast of struct: 
    vec_t v = vec_new(nel, sizeof(SOGrid_Node *));
    SOGrid_NodeRefVec r;
    r.nel = v.nel; r.el = (SOGrid_Node **)v.el;
    return r;
  }
*/

----------------------------------------------------------------------
    /* ** DEBUG ** */
    //    if(evalv[0] > 1)printf(" DEBUG!-> X: %16g Y: %16g \n", point[X], point[Y]);
    //    if(evalv[1] > 1)printf(" DEBUG!-> X: %16g Y: %16g \n", point[X], point[Y]);
    //    if(evalv[2] > 1)printf(" DEBUG!-> X: %16g Y: %16g \n", point[X], point[Y]);
    //    if(evalv[3] > 1)printf(" DEBUG!-> X: %16g Y: %16g \n", point[X], point[Y]);
    //    if(evalv[4] > 1)printf(" DEBUG!-> X: %16g Y: %16g \n", point[X], point[Y]);
/*     if((min[X] < 0.35)||(max[X] > 0.65))return; */
/*     if((min[Y] < 0.35)||(max[Y] > 0.65))return; */
/*     if(rank == 9)if(k != 754 && k != 752)return;  */
/*        printf(" DBG -> Index: %d Rank: %d \n", k, rank); */

void SO2DPlot_Begin_Figure
  (
    FILE *plot_file,
    double *min,
    double *max
  )
/*
  Prepares an Encapsulated Postscript file for drawing a dyadic SOGrid. */
{
  double mrg = 18.0;           // Margin (pt)
  double totd = 600.0;         // Diagonal (pt)

  double dx = max[X] - min[X];
  double dy = max[Y] - min[Y];
  
  double dd = hypot(dx,dy);    // Diagonal (user units).
  double umrg = mrg*(dd/totd); // Margin (user units)

  double toth = dx/dd*600.0 + 2*mrg;  // Total plot width (pt)
  double totv = dy/dd*600.0 + 2*mrg;  // Total plot height (pt)
 
  eps_begin_figure(plot_file, 
     min[X]-umrg, max[X]+umrg, min[Y]-umrg, max[Y]+umrg, 
     0.0, toth, 0.0, totv,
     1,1
  );
}


void SO2DPlot_End_Figure(FILE *plot_file)
/*
  Terminates an Encapsulated Postscript drawing of a dyadic SOGrid. */
{
  eps_end_figure(plot_file);
}



void SOProcFunction_AddMth(T *f, double a, T *h)
  { if (strcmp(f->type, h->type) != 0)
      { fprintf (stderr, "copy: type mismatch: \"%s\", \"%s\"\n", 
          f->type, h->type);
        affirm(FALSE, "aborted");
      }
    f->d->scale += h->d->scale;
  }

void SOProcFunction_ScaleMth(T *f, double a)
  { f->d->scale *= a;
  }

void SOProcFunction_MapleMth(T *f, FILE *wr)
  { affirm(FALSE , "maple not implemented yet");
  }

    fprintf(wr, "scale = %.16g\n",  f->d->scale);
    
      m->fn.maple = (MapleMth *)&SOProcFunction_MapleMth;
      m->fn.scale = (ScaleMth *)&SOProcFunction_ScaleMth;
      m->fn.add = (AddMth *)&SOProcFunction_AddMth;
  f->d->scale = 1.0;
    scale = nget_double(rd, "scale"); fget_eol(rd);
    f->d->scale = scale;
    
typedef void ScaleMth(OBJ *f, double a);       /* Scales function {f} by {a}. */
typedef void AddMth(OBJ *f, double a, OBJ *h); /* {f += a*h} (if compatible). */

ScaleMth *scale;   /* Scales the function {f} by {a}. */
AddMth *add;       /* Adds {a*h} to {f} (if compatible types). */



		//         printf(" Index: %d, lineval: %f, minp: %f, step_range: %f \n", k, lineval, minp, step_range);

		//         printf(" Px: %f, Py: %f \n",tri_vec.el[i].P[X], tri_vec.el[i].P[Y]);
		//         printf(" Qx: %f, Qy: %f \n",tri_vec.el[i].Q[X], tri_vec.el[i].Q[Y]);
		//         printf(" Rx: %f, Ry: %f \n",tri_vec.el[i].R[X], tri_vec.el[i].R[Y]);

		//         printf(" tri: %d, Pf: %f, Qf: %f, Rf: %f \n\n", i, trival.P, trival.Q, trival.R );


                //**************************** Debug triangles
	        //double tx[3], ty[3];
                //tx[0]=tri_vec.el[i].P[X]; tx[1]=tri_vec.el[i].Q[X]; tx[2]=tri_vec.el[i].R[X];
	        //ty[0]=tri_vec.el[i].P[Y]; ty[1]=tri_vec.el[i].Q[Y]; ty[2]=tri_vec.el[i].R[Y];
                //if(i == 0)eps_fill_polygon(plot_file, tx, ty, 3, 0.0, 0.0, 1.0);
                //if(i == 2)eps_fill_polygon(plot_file, tx, ty, 3, 0.0, 1.0, 0.0);
                //if(i == 4)eps_fill_polygon(plot_file, tx, ty, 3, 1.0, 0.0, 0.0);
                //if(i == 7)eps_fill_polygon(plot_file, tx, ty, 3, 0.0, 0.0, 0.0);
                //****************************************

	        //**************************** Debug cells

	        //if(k == 8)eps_fill_rectangle(plot_file,min[X],max[X],min[Y],max[Y],0.0,0.0,1.0);
	        //if(k == 10)eps_fill_rectangle(plot_file,min[X],max[X],min[Y],max[Y],1.0,0.0,0.0);
	        //if(k == 15)eps_fill_rectangle(plot_file,min[X],max[X],min[Y],max[Y],0.0,1.0,0.0);
	        //if(k == 13)eps_fill_rectangle(plot_file,min[X],max[X],min[Y],max[Y],0.0,0.0,0.0); 
                //****************************************


          

    double temp, tri_x[3], tri_y[3], sqr_x[4], sqr_y[4];
    double SPR, SRP, PRx, PRy, SQR, SRQ, QRx, QRy;
    int c_ind;

  /*   printf(" DBG ### >> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
  /*   printf(" DBG ### >> stp_rng: %f, color_ind: %d \n", fStep, color_ind); */

    if(!((Pf < 0.0)&&(Rf > 0.0)))
      { 
        if(line)
          {  if((Pf == 0.0)&&(Qf == 0.0))
               eps_draw_segment(psf, Px, Py, Qx, Qy); 

             if((Qf == 0.0)&&(Rf == 0.0))
               eps_draw_segment(psf, Qx, Qy, Rx, Ry); 

             if((Pf == 0.0)&&(Rf == 0.0))
               eps_draw_segment(psf, Px, Py, Rx, Ry); 
          }  
        return;
      }

    /* Now Pf must be negative and Rf must be positive */

    /* Swap q,r so that sign(Qf) != sign(Pf): */

    if( Qf < 0.0 )
      {
        temp = Rx; Rx = Px; Px = temp;
        temp = Ry; Ry = Py; Py = temp;
        temp = Rf; Rf = Pf; Pf = temp; 
      }

    affirm( Rf != 0.0, "Rf == 0.0" );
    affirm( Pf != 0.0, "Pf == 0.0" );
    affirm((Rf > 0.0) != (Pf > 0.0), "Rf sign == Pf sign" );
    affirm((Qf == 0.0)||(((Qf > 0.0) == (Rf > 0.0))&&((Qf > 0.0) != (Pf > 0.0))), "Qf sign == Pf sign" );

    if( Qf != 0.0 )
      {

        SQR = (double)Pf / (Pf - Rf);
        SRQ = (double)Rf / (Rf - Pf);

        QRx = (double)Rx * SQR + (double)Px * SRQ;
        QRy = (double)Ry * SQR + (double)Py * SRQ;

        SPR = (double)Pf / (Pf - Qf);
        SRP = (double)Qf / (Qf - Pf);

        PRx = (double)Qx * SPR + (double)Px * SRP;
        PRy = (double)Qy * SPR + (double)Py * SRP;

        if(line){eps_draw_segment(psf, PRx, PRy, QRx, QRy); return;}

        int k;
        for(k=-1; k < 2; k=k+2)
          { 
            if(k * Pf > 0)
              {       
                if(Pf > 0.0) c_ind = color_ind + 1; else c_ind = color_ind; 
  /*               if(c_ind < 0 || c_ind >= steps+2) */
  /*                 {printf(" DBG!! c_ind: %d, color_ind: %d \n", c_ind, color_ind);  */
  /*                  c_ind = steps+1;}            */

                if((fabs(Pf) <= fStep)||((color_ind == 0)&&(Pf < 0)))
                  {
                    tri_x[0] = PRx; tri_x[1] = QRx; tri_x[2] = Px;
                    tri_y[0] = PRy; tri_y[1] = QRy; tri_y[2] = Py;
                    eps_fill_polygon(psf, tri_x, tri_y, 3, 
                        colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /* 		  printf(" DBG 0>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
                   }
                else
                  { if(Pf < 0)
                      { sqr_x[0] = PRx; sqr_x[1] = QRx; 
                        sqr_x[2] = plot_p0[0]; sqr_x[3] = plot_p1[0];
                        sqr_y[0] = PRy; sqr_y[1] = QRy; 
                        sqr_y[2] = plot_p0[1]; sqr_y[3] = plot_p1[1];
                        eps_fill_polygon(psf, sqr_x, sqr_y, 4, 
                           colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]); 
  /*                       printf(" DBG 1>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */

                        if((fabs(Rf) > fStep)&&(fabs(Qf) < fStep))
                          { tri_x[0] = PRx; tri_x[1] = plot_p1[0]; tri_x[2] = Qx;
                            tri_y[0] = PRy; tri_y[1] = plot_p1[1]; tri_y[2] = Qy;  
                            eps_fill_polygon(psf, tri_x, tri_y, 3, 
                               colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*                           printf(" DBG 2>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
                          }
                      }    
                    else
                      { plot_p0[0] = QRx; plot_p1[0] = PRx;
                        plot_p0[1] = QRy; plot_p1[1] = PRy;
  /*                    printf(" DBG 3>> Qf: %.16g, Rf: %.16g, Pf: %.16g \n", Qf, Rf, Pf); */
                      } 
                  }            
              }   

            if(k * Rf > 0) 
              {   
                if(Rf > 0.0) c_ind = color_ind + 1; else c_ind = color_ind;
  /*               if(c_ind < 0 || c_ind >= steps+2) */
  /*                 {printf(" DBG!! c_ind: %d, color_ind: %d \n", c_ind, color_ind); */
  /*                  c_ind = steps+1;} */

                if(((fabs(Qf) <= fStep)&&(fabs(Rf) <= fStep))||((color_ind == 0)&&(Rf < 0)))
                  { 
                    sqr_x[0] = Qx; sqr_x[1] = Rx; sqr_x[2] = QRx; sqr_x[3] = PRx;
                    sqr_y[0] = Qy; sqr_y[1] = Ry; sqr_y[2] = QRy; sqr_y[3] = PRy;
                    eps_fill_polygon(psf, sqr_x, sqr_y, 4, 
                        colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*                   printf(" DBG 4>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
                  }  
                else
                  { if(Rf < 0)
                      { sqr_x[0] = PRx; sqr_x[1] = QRx; 
                        sqr_x[2] = plot_p0[0]; sqr_x[3] = plot_p1[0];
                        sqr_y[0] = PRy; sqr_y[1] = QRy; 
                        sqr_y[2] = plot_p0[1]; sqr_y[3] = plot_p1[1];
                        eps_fill_polygon(psf, sqr_x, sqr_y, 4, 
                            colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*                       printf(" DBG 5>> Qf: %.16g, Rf: %.16g, Pf: %.16g \n", Qf, Rf, Pf); */

                        if((fabs(Rf) > fStep)&&(fabs(Qf) < fStep))
                          { tri_x[0] = PRx; tri_x[1] = plot_p1[0]; tri_x[2] = Qx;
                            tri_y[0] = PRy; tri_y[1] = plot_p1[1]; tri_y[2] = Qy;  
                            eps_fill_polygon(psf, tri_x, tri_y, 3, 
                               colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*                           printf(" DBG 6>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
                          } 
                      }        
                    else
                      { plot_p0[0] = QRx; plot_p1[0] = PRx;                
                        plot_p0[1] = QRy; plot_p1[1] = PRy;   
  /*                       printf(" DBG 7>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
                      } 
                  }   
              }
          }     
      }
    else
      { 
        SQR = (double)Pf / (Pf - Rf);
        SRQ = (double)Rf / (Rf - Pf);
        QRx = (double)Rx * SQR + (double)Px * SRQ;
        QRy = (double)Ry * SQR + (double)Py * SRQ;

        if(line){eps_draw_segment(psf, Qx, Qy, QRx, QRy); return;} 

        if(Pf > 0.0) c_ind = color_ind + 1; else c_ind = color_ind;
  /*       if(c_ind < 0 || c_ind >= steps+2) */
  /*         {printf(" DBG!! c_ind: %d, color_ind: %d \n", c_ind, color_ind); */
  /*          c_ind = steps+1;} */

        if((fabs(Pf) <= fStep)||((color_ind == 0)&&(Pf < 0)))
          {
            tri_x[0] = Qx; tri_x[1] = Px; tri_x[2] = QRx;
            tri_y[0] = Qy; tri_y[1] = Py; tri_y[2] = QRy;
            eps_fill_polygon(psf, tri_x, tri_y, 3, 
               colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*           printf(" DBG 8>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
          }
        else
          { if(Pf < 0)
              { sqr_x[0] = Qx; sqr_x[1] = QRx; 
                sqr_x[2] = plot_p0[0]; sqr_x[3] = plot_p1[0];
                sqr_y[0] = Qy; sqr_y[1] = QRy; 
                sqr_y[2] = plot_p0[1]; sqr_y[3] = plot_p1[1];
                eps_fill_polygon(psf, sqr_x, sqr_y, 4, 
                   colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*               printf(" DBG 9>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */

                if((fabs(Rf) > fStep)&&(fabs(Qf) < fStep))
                  { tri_x[0] = PRx; tri_x[1] = plot_p1[0]; tri_x[2] = Qx;
                    tri_y[0] = PRy; tri_y[1] = plot_p1[1]; tri_y[2] = Qy;  
                    eps_fill_polygon(psf, tri_x, tri_y, 3, 
                       colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*                   printf(" DBG 10>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
                  }
              }
            else
              { plot_p0[0] = QRx; plot_p1[0] = Qx;
                plot_p0[1] = QRy; plot_p1[1] = Qy; 
  /*               printf(" DBG 11>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
              }
          }

        if(Rf > 0.0) c_ind = color_ind + 1; else c_ind = color_ind;
  /*       if(c_ind < 0 || c_ind >= steps+2) */
  /*         {printf(" DBG!! c_ind: %d, color_ind: %d \n", c_ind, color_ind); */
  /*          c_ind = steps+1;} */

        if(fabs(Rf) <= fStep)
          {
            tri_x[0] = Qx; tri_x[1] = Rx; tri_x[2] = QRx;
            tri_y[0] = Qy; tri_y[1] = Ry; tri_y[2] = QRy;
            eps_fill_polygon(psf, tri_x, tri_y, 3, 
               colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*           printf(" DBG 12>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
          }
        else
          { if(Rf < 0)
              { sqr_x[0] = Qx; sqr_x[1] = QRx; 
                sqr_x[2] = plot_p0[0]; sqr_x[3] = plot_p1[0];
                sqr_y[2] = plot_p0[1]; sqr_y[3] = plot_p1[1];
                eps_fill_polygon(psf, sqr_x, sqr_y, 4, 
                   colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*               printf(" DBG 13>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */

                if((fabs(Rf) > fStep)&&(fabs(Qf) < fStep))
                  { tri_x[0] = PRx; tri_x[1] = plot_p1[0]; tri_x[2] = Qx;
                    tri_y[0] = PRy; tri_y[1] = plot_p1[1]; tri_y[2] = Qy;  
                    eps_fill_polygon(psf, tri_x, tri_y, 3, 
                       colors[c_ind][0], colors[c_ind][1], colors[c_ind][2]);
  /*                   printf(" DBG 14>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
                  }
              }
            else
              { plot_p0[0] = QRx; plot_p1[0] = Qx;
                plot_p0[1] = QRy; plot_p1[1] = Qy;
  /*               printf(" DBG 15>> Qf: %f, Rf: %f, Pf: %f \n", Qf, Rf, Pf); */
              }
          }

      }
  }


#define Rank_Inc 4
#define Max_Rank 10
#define Tri_Lines 15
double Iso_min, Iso_max; /* minimun and maximun values of plotted function */

SO2DPlot_Triangle_vec SO2DPlot_Cell_Triangles
  ( double *min,
    double *max
  );
  /* Returns a vector of 8 triangles (SO2DPlot_Triangle_vec) 
     for a given cell, each one with the center and 2 cell 
     vertices points */

void SO2DPlot_Triangle_Plot
  (
    FILE* psf, 	
    SO2DPlot_Triangle* tri
  );


void SO2DPlot_Triangle_0_Plot
  (
    FILE* psf, 		
    double Px, double Py, double Pf,
    double Qx, double Qy, double Qf,
    double Rx, double Ry, double Rf,
    int color_range,
    int steps,                        /* Number of isolines between {fPlotMin} and {fPlotMax} */
    bool_t line,
    double fStep,
    double *plot_p0, 
    double *plot_p1 
  );
  /*
    Given a triangle "P,Q,R" on the plane, and function values at its
    vertices, interpolates a linear function through that data, and
    plots the set of points where that function is zero. */

  
void SO2DPlot_Plot_Subtree(
  FILE *psf, 
  SOGrid_Node *p, 
  double *min,
  double *max,
  int long_axis,
  int depth,
  int max_depth
);
   
    if(fPlotMin < 0)fsat = fPlotMin;
    fvalue = fPlotMin; 

    while((fvalue < 0)&&(ind < (steps+1)))
      { 
        t = (double)fvalue / fsat;

/*         printf(" DBG: IND: %d, T: %.16g, fvalue: %.16g, fsat: %.16g, \n", ind, t, fvalue, fsat); */

        bandColor[ind][0] = fabs(1 - t);  /* red */
        bandColor[ind][1] = fabs(1 - t);  /* green */
        bandColor[ind][2] = 1;            /* blue */

/*         printf(" -> c0: %.16g, c1: %.16g, c2: %.16g \n\n", */
/*            bandColor[ind][0],bandColor[ind][1],bandColor[ind][2]); */

        ind++; fvalue += fStep;
      }

    if(fPlotMax >= 0)
      { fsat = fPlotMax;
        bandColor[ind][0] = 1; /* red */
        bandColor[ind][1] = 1; /* green */
        bandColor[ind][2] = 1; /* blue */

/*      printf("Zero-> DBG: IND: %d, T: %.16g, fvalue: %.16g, fsat: %.16g, \n", ind, t, fvalue, fsat); */
/*      printf("Zero -> c0: %.16g, c1: %.16g, c2: %.16g \n\n", */
/*              bandColor[ind][0],bandColor[ind][1],bandColor[ind][2]); */

        ind++; fvalue += fStep;
      }

    while((fvalue > 0)&&(ind < (steps+1)))
      { 
        t = (double)fvalue / fsat;

/*         printf(" DBG: IND: %d, T: %.16g, fvalue: %.16g, fsat: %.16g, \n", ind, t, fvalue, fsat); */

        bandColor[ind][0] = 1;           /* red */
        bandColor[ind][1] = fabs(1 - t); /* green */
        bandColor[ind][2] = fabs(1 - t); /* blue */

/*         printf(" -> c0: %.16g, c1: %.16g, c2: %.16g \n\n", */
/*         bandColor[ind][0],bandColor[ind][1],bandColor[ind][2]); */

        ind++; fvalue += fStep;
      }

void SO2DPlot_Triangle_Plot
  (
    FILE* psf, 	
    SO2DPlot_Triangle* tri
  )
  {
    eps_set_pen(psf, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);

    eps_draw_segment(psf, tri->P[X], tri->P[Y], tri->Q[X], tri->Q[Y]);
    eps_draw_segment(psf, tri->P[X], tri->P[Y], tri->R[X], tri->R[Y]);
    eps_draw_segment(psf, tri->R[X], tri->R[Y], tri->Q[X], tri->Q[Y]);

    eps_set_pen(psf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
   
  }
========================================

    /* Plot coordinates: */
    double minPlot[PDIM], maxPlot[PDIM]; /* Client plot ranges for {C}. */

    double minp, maxp, lineval;
    double plot_p0[PDIM], plot_p1[PDIM]; 
    SO2DPlot_Triangle_Val trival;   
    int color_ind; 


====================
            
            
        SO2DPlot_Triangle_vec tri_vec, tri_vecPlot; 
        tri_vecPlot =  SO2DPlot_Cell_Triangles(minPlot, maxPlot);
        tri_vec =  SO2DPlot_Cell_Triangles(minRel, maxRel);

        for(i = 0; i < tri_vec.nel; i++)
          { 
            /* Gets triangle values (for P Q R) */
            f->m->eval(f, tri_vec.el[i].P, &trival.P); //if(trival.P < 0)trival.P = 0;
            f->m->eval(f, tri_vec.el[i].Q, &trival.Q); //if(trival.Q < 0)trival.Q = 0; 
            f->m->eval(f, tri_vec.el[i].R, &trival.R); //if(trival.R < 0)trival.R = 0; 
 
/*    	    if(k > 650) continue; */
/*          if(rank != 9)continue; */
/* 	    if(i != 0)continue; */

            minp = trival.P; 
            if(trival.Q < minp)minp = trival.Q; if(trival.R < minp)minp = trival.R;        
            maxp = trival.P; 
            if(trival.Q > maxp)maxp = trival.Q; if(trival.R > maxp)maxp = trival.R;

            color_ind = (int)((minp - (fPlotMin + (double)fStep/2)) / fStep); 
            if(color_ind < 0)color_ind = 0; if(color_ind >= steps)color_ind = steps-1;
            lineval = (double)(color_ind * fStep) + (fPlotMin + (double)fStep/2);

/*             printf(" DBG## minp: %.16g, fPlotMin: %.16g, fStep: %.16g \n", */
/*                      minp, fPlotMin, fStep); */
            
/*             printf("DBG lineval: %.16g, color_ind: %d \n", lineval, color_ind); */
/*             printf("DBG Qf: %.16g, Rf: %.16g, Pf: %.16g \n\n", trival.P, trival.Q, trival.R); */
            //printf(" DBG# steps: %d \n", steps);

            if((minp >= lineval)&&(maxp <= fStep + lineval) && !line)
              { double tri_x[3], tri_y[3];

/* 		if(color_ind < 0 || color_ind >= steps+2) */
/*              { printf(" DBG!! c_ind: %d, color_ind: %d steps: %d \n", c_ind, color_ind, steps); */
/*                c_ind = steps+1; } */

                tri_x[0] = tri_vec_corr.el[i].P[X];
                tri_x[1] = tri_vec_corr.el[i].Q[X];
                tri_x[2] = tri_vec_corr.el[i].R[X];
                tri_y[0] = tri_vec_corr.el[i].P[Y];
                tri_y[1] = tri_vec_corr.el[i].Q[Y];
                tri_y[2] = tri_vec_corr.el[i].R[Y];

/*                 printf(" DBG## color_ind: %d, steps: %d \n\n", color_ind, steps); */

                eps_fill_polygon(psf, tri_x, tri_y, 3,
		    bandColor[color_ind+1][0], bandColor[color_ind+1][1], bandColor[color_ind+1][2]);
             
/* 		printf(" DBG!! c0: %.16g, c1: %.16g, c2: %.16g \n", */
/* 		      bandColor[color_ind+1][0], bandColor[color_ind+1][1], bandColor[color_ind+1][2]); */

                continue;
              }
   
	    if(maxp > fPlotMax)maxp = fPlotMax; j = 0;

            while((lineval <= maxp)&&(j < Tri_Lines))
              {
                /* printf("\n\n"); */
/*                 printf("DBG lineval: %.16g, j: %d, color_ind: %d \n", lineval, j, color_ind); */
/*                 printf("DBG fPlotMin: %.16g, fPlotMax: %.16g step_rng: %.16g \n", fPlotMin, fPlotMax, fStep); */

                eps_set_pen(psf, 1.0, 0.0, 0.0, 0.2, 0.0, 0.0);
                SO2DPlot_Triangle_0_Plot /* Plots zero line inside triangle */
                  (
                    psf,
                    tri_vec_corr.el[i].P[X], tri_vec_corr.el[i].P[Y],trival.P - lineval,
                    tri_vec_corr.el[i].Q[X], tri_vec_corr.el[i].Q[Y],trival.Q - lineval,
                    tri_vec_corr.el[i].R[X], tri_vec_corr.el[i].R[Y],trival.R - lineval,
                    color_ind, steps, line, fStep, plot_p0, plot_p1
                    
                  );
                eps_set_pen(psf, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

                color_ind++; j++;
                /* lineval represents the value of the isoline to be plotted */
                lineval += fStep;
              }  
            //SO2DPlot_Triangle_Plot(psf, &tri_vec_corr.el[i]); /* Plots triangles */
          }   
      }

  }


/* Returns a vector of 8 triangles (SO2DPlot_Triangle_vec) 
   for a given cell, each one with the center and 2 cell 
   vertices points */

SO2DPlot_Triangle_vec SO2DPlot_Cell_Triangles
  ( double *min,
    double *max
  )
  {
    double center[PDIM];
    SO2DPlot_Triangle_vec tri_vec = SO2DPlot_Triangle_vec_new(8);

    /* Gets the center point coordinates */
    center[X] = (min[X] + max[X]) * 0.5;
    center[Y] = (min[Y] + max[Y]) * 0.5;
    
    /* Sets the coordinates for triangle 1 */ 
    tri_vec.el[0].P[X] = min[X]; 
    tri_vec.el[0].P[Y] = min[Y];

    tri_vec.el[0].Q[X] = center[X];
    tri_vec.el[0].Q[Y] = min[Y];

    tri_vec.el[0].R[X] = center[X];
    tri_vec.el[0].R[Y] = center[Y];

    /* Sets the coordinates for triangle 2 */
    tri_vec.el[1].P[X] = center[X];
    tri_vec.el[1].P[Y] = min[Y];

    tri_vec.el[1].Q[X] = max[X];
    tri_vec.el[1].Q[Y] = min[Y];

    tri_vec.el[1].R[X] = center[X];
    tri_vec.el[1].R[Y] = center[Y];

    /* Sets the coordinates for triangle 3 */
    tri_vec.el[2].P[X] = max[X];
    tri_vec.el[2].P[Y] = min[Y];

    tri_vec.el[2].Q[X] = max[X];
    tri_vec.el[2].Q[Y] = center[Y];

    tri_vec.el[2].R[X] = center[X];
    tri_vec.el[2].R[Y] = center[Y];

    /* Sets the coordinates for triangle 4 */
    tri_vec.el[3].P[X] = max[X];
    tri_vec.el[3].P[Y] = center[Y];

    tri_vec.el[3].Q[X] = max[X];
    tri_vec.el[3].Q[Y] = max[Y];

    tri_vec.el[3].R[X] = center[X];
    tri_vec.el[3].R[Y] = center[Y];

    /* Sets the coordinates for triangle 5 */
    tri_vec.el[4].P[X] = max[X];
    tri_vec.el[4].P[Y] = max[Y];

    tri_vec.el[4].Q[X] = center[X];
    tri_vec.el[4].Q[Y] = max[Y];

    tri_vec.el[4].R[X] = center[X];
    tri_vec.el[4].R[Y] = center[Y];

    /* Sets the coordinates for triangle 6 */
    tri_vec.el[5].P[X] = center[X];
    tri_vec.el[5].P[Y] = max[Y];

    tri_vec.el[5].Q[X] = min[X];
    tri_vec.el[5].Q[Y] = max[Y];

    tri_vec.el[5].R[X] = center[X];
    tri_vec.el[5].R[Y] = center[Y];
 
    /* Sets the coordinates for triangle 7 */
    tri_vec.el[6].P[X] = min[X];
    tri_vec.el[6].P[Y] = max[Y];

    tri_vec.el[6].Q[X] = min[X];
    tri_vec.el[6].Q[Y] = center[Y];

    tri_vec.el[6].R[X] = center[X];
    tri_vec.el[6].R[Y] = center[Y];

    /* Sets the coordinates for triangle 8 */  
    tri_vec.el[7].P[X] = min[X];
    tri_vec.el[7].P[Y] = center[Y];

    tri_vec.el[7].Q[X] = min[X];
    tri_vec.el[7].Q[Y] = min[Y];

    tri_vec.el[7].R[X] = center[X];
    tri_vec.el[7].R[Y] = center[Y];

    return tri_vec;
}
  
void SO2DPlot_Plot_Tree(
  FILE *psf, 
  SOGrid_Tree *tree, 
  double *min,
  double *max,
  int max_depth
)
{
  affirm(tree->d == 2, "wrong tree dimension");
  eps_set_pen(psf, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0);
  eps_draw_rectangle(psf, min[X], max[X], min[Y], max[Y]);
  SO2DPlot_Plot_Subtree(psf, tree->root, min, max, X, 1, max_depth);
}
  
void SO2DPlot_Plot_Subtree(
  FILE *psf, 
  SOGrid_Node *p, 
  double *min,
  double *max,
  int long_axis,
  int depth,
  int max_depth
)
{
  double min_aux[PDIM], max_aux[PDIM];

  min_aux[X] = min[X];
  max_aux[X] = max[X];
  min_aux[Y] = min[Y];
  max_aux[Y] = max[Y];

  if ((p != NULL) && ((p->c[0] != NULL) || (p->c[1] != NULL)) && depth < max_depth)
    { 
      double xmid = (max[X] + min[X])/2.0;
      double ymid = (max[Y] + min[Y])/2.0;

      eps_set_pen(psf, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0);
      //      eps_fill_dot(psf, xmid, ymid, 0.5, 0, 0, 0);

      if (long_axis == X)
      { 
          eps_draw_segment(psf, xmid, min[Y], xmid, max[Y]);
	  max_aux[X] = xmid;
          SO2DPlot_Plot_Subtree(psf, p->c[0], min_aux, max_aux, Y, depth + 1, max_depth);

	  max_aux[X] = max[X];
	  min_aux[X] = xmid;
          SO2DPlot_Plot_Subtree(psf, p->c[1], min_aux, max_aux, Y, depth + 1, max_depth);
      }
      else
      { 
          eps_draw_segment(psf, min[X], ymid, max[X], ymid);
	  max_aux[Y] = ymid;
          SO2DPlot_Plot_Subtree(psf, p->c[0], min_aux, max_aux, X, depth + 1, max_depth);

	  max_aux[Y] = max[Y];
	  min_aux[Y] = ymid;

          SO2DPlot_Plot_Subtree(psf, p->c[1], min_aux, max_aux, X, depth + 1, max_depth);
      }
    }
}
