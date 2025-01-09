/* See {tmaze_figs_bb.h}. */
/* Last edited on 2024-12-21 11:31:36 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <bool.h>
#include <jsstring.h>
#include <dgraph.h>
#include <jsrandom.h>

#include <tmaze.h>
#include <tmaze_bb.h>
#include <tmaze_bb_graph.h>
#include <tmaze_figs_bb.h>

#define tmaze_bb_CURVED TRUE

void tmaze_bb_plot_pattern_tileset_curved(epswr_figure_t *eps, tmaze_pattern_t *M, int x, int y, bool_t fill, bool_t draw);
  /* Plots the tileset in column {x} and  row {y} of maze pattern {M}, with default style. */

void tmaze_bb_plot_pattern_tileset_curved(epswr_figure_t *eps, tmaze_pattern_t *M, int x, int y, bool_t fill, bool_t draw)
  { tmaze_bb_plot_pattern_tileset(eps, M, x, y, tmaze_bb_CURVED, fill, draw); }

void tmaze_figs_plt_bb_tileset(char *dir, char *name, tmaze_tileset_t tset)
  {
    tmaze_pattern_t M = tmaze_make_pattern(1, 1, FALSE, TRUE, 1);
    tmaze_fill_pattern(&M, tset);
    
    tmaze_plot_tileset_proc_t *proc = &tmaze_bb_plot_pattern_tileset_curved;
    tmaze_plot_pattern(txtcat("out/",dir), name, &M, TRUE, TRUE, proc, FALSE, 12.0, 1.0);
  }
  
void tmaze_figs_plt_bb_sample(char *dir, char *name, int nx, int ny, int seed, bool_t showbig)
  {

    auto void plot_tile(epswr_figure_t *eps, tmaze_t *M, int x, int y, bool_t fill, bool_t draw);
      /* Plots the tile in column {x} and  row {y} of maze {M}, with default style. */

    auto frgb_t *pick_color(int e);
      /* Picks a color for edge {e}. The edge numbering is as in
        {tmaze_bb_graph_get_edge_components}. */

    /* Buil the maze and choose the tile size: */
    tmaze_t M;
    double tsize;
    bool_t debug;
    if ((nx == 6) && (ny == 5))
      { M = tmaze_bb_from_string
          ( ">>>>>>;"
            ">><>>>;"
            "><><>>;"
            ">>>><>;"
            ">>><>>;",
            TRUE
          );
        tsize = 8.0;
        debug = TRUE;
      }
    else if ((nx == 12) && (ny == 7))
      { M = tmaze_bb_from_string
          ( ">>><<<>><<<<;"
            "><<>><<<><<<;"
            "<<><<><><<>>;"
            ">><><<<>>>>>;"
            ">><>>>><>><<;"
            ">>><<><><><>;"
            "<<<<<<><<<<<;",
            TRUE
          );
        tsize = 8.0;
        debug = FALSE;
      }
    else
      { M = tmaze_bb_make(nx, ny, TRUE); 
        tmaze_bb_random(&M, seed);
        tsize = 4.233;
        debug = FALSE;
      }
      
    tmaze_cell_count_t nt = M.nt;
    
    /* These variables are set only if {showbig} has been requested: */
    dgraph_vertex_count_t nv = 0;      /* Number of graph vertices. */
    dgraph_edge_count_t ne = 0;        /* Number of graph edges. */
    tmaze_comp_count_t nc = 0;         /* Number of connected components. */
    tmaze_size_t *sz = NULL;           /* V-count of component of each vertex {v}. */
    dgraph_vertex_index_t *rt = NULL;  /* Root vertex of each edge {e}. */
    dgraph_vertex_index_t rmax1 = nv;  /* Root of largest component. */
    dgraph_vertex_index_t rmax2 = nv;  /* Root of second largest component. */
    if (showbig)
      { /* Get the graph in order to find the big components:*/
        dgraph_t G = tmaze_bb_graph_make(&M);
        assert(G.rows == G.cols);
        nv = G.rows;

        /* Get the size of the component that contains each edge: */
        tmaze_bb_graph_get_edge_components(&G, &M, &ne, &nc, &sz, &rt);
        
        /* Find largest and secon-largest components {rmax1,rmax2}: */
        rmax1 = nv;  /* Root of largest component. */
        rmax2 = nv;  /* Root of second largest component. */
        int e;
        for (e = 0; e < ne; e++)
          { /* Get root of edge {e} and its size: */
            dgraph_vertex_index_t re = rt[e]; assert(re < nv);
            tmaze_size_t sze = sz[re];
            /* Update {rmax1,rmax2}: */
            if ((rmax1 >= nv) || (sze > sz[rmax1]))
              { rmax2 = rmax1; rmax1 = re; }
            else if ((re != rmax1) && ((rmax2 >= nv) || (sze > sz[rmax2])))
              { rmax2 = re; }
          }
        assert(rmax1 < nv);
        assert(rmax2 < nv);
        fprintf(stderr, "%9d tiles\n", nt);
        fprintf(stderr, "%9d vertices\n", nv);
        fprintf(stderr, "%9d edges\n", ne);
        fprintf(stderr, "%9d components\n", nc);
        fprintf(stderr, "component with root %3d has %9d vertices\n", rmax1, sz[rmax1]);
        fprintf(stderr, "component with root %3d has %9d vertices\n", rmax2, sz[rmax2]);

        /* Count and print components by size, for documentation: */
        tmaze_size_t ms;
        tmaze_comp_count_t *ct; /* Count cmponents by size. */
        tmaze_bb_graph_count_components_by_size(&G, &M, &ms, &ct);
        assert(ms <= nv);
        fprintf(stderr, "  %9s %9s\n", "size", "comps");
        fprintf(stderr, "  %9s %9s\n", "---------", "---------");
        tmaze_size_t s;
        for (s = 0; s <= ms; s++)
          { if (ct[s] != 0) { fprintf(stderr, "  %9d %9d\n", s, ct[s]); } }
        fprintf(stderr, "  %9s %9s\n", "---------", "---------");

        /* Cleanup: */
        dgraph_trim(&G, 0);
        free(ct);
      }
    
    /* Colors for road filling: */
    frgb_t rgb_lite = (frgb_t){{ 1.0000f, 0.9500f, 0.6000f }}; /* Light RGB color. */
    frgb_t rgb_mean = (frgb_t){{ 0.5000f, 0.8500f, 1.0000f }}; /* Medium RGB color. */
    frgb_t rgb_dark = (frgb_t){{ 1.0000f, 0.5500f, 0.2000f }}; /* Strong RGB color. */
    
    tmaze_plot_tile_proc_t *proc = &plot_tile;
    tmaze_plot_maze(txtcat("out/",dir), name, &M, TRUE, TRUE, proc, FALSE, tsize, 1.0);
    
    /* Cleanup: */
    free(sz);
    return;

    /* INTERNAL PROCS */
    
    frgb_t *pick_color(int e)
      {
        if (rt == NULL) { return NULL; }
        dgraph_vertex_index_t re = rt[e]; assert(re < nv);
        if (re == rmax1)
          { return &rgb_dark; }
        else if (re == rmax2)
          { return &rgb_mean; }
        else  
          { return &rgb_lite; }
      }
    
    void plot_tile(epswr_figure_t *eps, tmaze_t *M, int x, int y, bool_t fill, bool_t draw)
      { 
        frgb_t *rgb_S = NULL; /* Color for the South road. */
        frgb_t *rgb_N = NULL; /* Color for the North road. */
        if (showbig)
          { assert(rt != NULL);
            assert(sz != NULL);
            /* Find the components of the South and North edges: */
            tmaze_cell_index_t k = tmaze_cell_index(x, y, M->nx, M->ny);
            int eS = 2*k; assert(eS < ne);
            int eN = 2*k + 1; assert(eN < ne);
            if (debug)
              { fprintf(stderr, "tile (%2d,%2d) = %4d", x, y, k);
                fprintf(stderr, "  rt[eS] = %5d  rt[eN] = %5d", rt[eS], rt[eN]);
                fprintf(stderr, "\n");
              }
            /* Choose the colors of the two edges: */
            rgb_S = pick_color(eS);
            rgb_N = pick_color(eN);
          }
        tmaze_bb_plot_maze_tile(eps, M, x, y, tmaze_bb_CURVED, rgb_S, rgb_N, fill, draw);
      }

  }
  
void tmaze_figs_plt_bb_conjuntive_pat(char *dir, char *name)
  {
    tmaze_pattern_t M;
    M = tmaze_bb_pattern_from_string
      ( "**>>*;"
        "<><><;"
        "*>>**;",
        FALSE, TRUE, 1
      );
    tmaze_plot_tileset_proc_t *proc = &tmaze_bb_plot_pattern_tileset_curved;
    tmaze_plot_pattern(txtcat("out/",dir), name, &M, TRUE, TRUE, proc, FALSE, 8.0, 1.0);
  }
  
void tmaze_figs_plt_bb_component_pattern(char *dir, char *prefix, int k, int which)
  {
    tmaze_pattern_t M;
    switch(k)
      {
      case  4: M = tmaze_figs_make_bb_component_pattern_04(which); break;
      case  8: M = tmaze_figs_make_bb_component_pattern_08(which); break;
      case 12: M = tmaze_figs_make_bb_component_pattern_12(which); break;
      case 16: M = tmaze_figs_make_bb_component_pattern_16(which); break;
      case 20: M = tmaze_figs_make_bb_component_pattern_20(which); break;
      default: 
        fprintf(stderr, "** %s not implemented for k = %d\n", __FUNCTION__, k);
        exit(1);
      }
      
    char *name = jsprintf("%s-%02d-%02d", prefix, k, which);
    tmaze_plot_tileset_proc_t *proc = &tmaze_bb_plot_pattern_tileset_curved;
    tmaze_plot_pattern(txtcat("out/",dir), name, &M, TRUE, TRUE, proc, FALSE, 8.0, 1.0);
    free(name);
  }
  
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_04(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_bb_pattern_from_string("<>;><;", FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }
  
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_08(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_bb_pattern_from_string("*<>;<<<;><*;", FALSE, TRUE, 1); break;
      case 1: M = tmaze_bb_pattern_from_string("<>*;>>>;*><;", FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }
  
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_12(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_bb_pattern_from_string("**<>;*<<<;<<<*;><**;", FALSE, TRUE, 1); break;
      case 1: M = tmaze_bb_pattern_from_string("*<>*;<<>>;><><;",      FALSE, TRUE, 1); break;
      case 2: M = tmaze_bb_pattern_from_string("*<>*;<<>>;>><<;*><*;", FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }

tmaze_pattern_t tmaze_figs_make_bb_component_pattern_16(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_bb_pattern_from_string("***<>;**<<<;*<<<*;<<<**;><***;", FALSE, TRUE, 1); break;
      case 1: M = tmaze_bb_pattern_from_string("*<>*;<<>>;><<<;*><*;",           FALSE, TRUE, 1); break;
      case 2: M = tmaze_bb_pattern_from_string("*<><>;<<><<;><><*;",             FALSE, TRUE, 1); break;
      case 3: M = tmaze_bb_pattern_from_string("<>**;>>>*;<<>>;><><;",           FALSE, TRUE, 1); break;
      case 4: M = tmaze_bb_pattern_from_string("<>**;>>>*;<<>>;>><<;*><*;",      FALSE, TRUE, 1); break;
      case 5: M = tmaze_bb_pattern_from_string("*<>**;<<>>*;>>*>>;*>><<;**><*;", FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d not implemented yet\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }

tmaze_pattern_t tmaze_figs_make_bb_component_pattern_20(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_bb_pattern_from_string("*<><>*;<<><>>;>>**<<;*>><<*;**><**;",   FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }
