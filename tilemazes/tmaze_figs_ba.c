/* See {tmaze_figs_ba.h}. */
/* Last edited on 2024-12-21 11:31:30 by stolfi */

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
#include <jsrandom.h>

#include <tmaze.h>
#include <tmaze_ba.h>
#include <tmaze_figs_ba.h>

void tmaze_figs_plt_ba_tileset(char *dir, char *name, tmaze_tileset_t tset)
  {
    tmaze_pattern_t M = tmaze_make_pattern(1, 1, FALSE, TRUE, 1);
    tmaze_fill_pattern(&M, tset);
    tmaze_plot_tileset_proc_t *proc = &tmaze_ba_plot_pattern_tileset;
    tmaze_plot_pattern(txtcat("out/",dir), name, &M, TRUE, TRUE, proc, FALSE, 12.0, 1.0);
  }
  
void tmaze_figs_plt_ba_sample(char *dir, char *name, int nx, int ny, int seed, bool_t showbig)
  {
    tmaze_t M;
    double tsize;
    if ((nx == 8) && (ny == 4))
      { M = tmaze_ba_from_string
          ( "WNSEWNSE;"
            "NSEWNSEW;"
            "SEWNSEWN;"
            "EWNSEWNS;",
            TRUE
          );
        tsize = 8.0;
      }
    else if ((nx == 12) && (ny == 7))
      { M = tmaze_ba_from_string
          ( "NSNEEENNWEEE;"
            "SEENSEEWNEEE;"
            "EENEENENEENS;"
            "SNESWWWSNSNN;"
            "SNENSNNWSSWE;"
            "NNSEESENESEN;"
            "EWEEEESWEWWW;",
            TRUE
          );
        tsize = 8.0;
      }
    else
      { M = tmaze_ba_make(nx, ny, TRUE);
        tmaze_ba_random(&M, seed);
        tsize = 4.233;
      }
    tmaze_plot_tile_proc_t *proc = &tmaze_ba_plot_maze_tile;
    tmaze_plot_maze(txtcat("out/",dir), name, &M, TRUE, TRUE, proc, FALSE, tsize, 1.0);
  }
  
void tmaze_figs_plt_ba_conjuntive_pat(char *dir, char *name)
  {
    tmaze_pattern_t M;
    M = tmaze_ba_pattern_from_string
      ( "**SS*;"
        "ENWSW;"
        "*SS**;",
        FALSE, TRUE, 1
      );
    tmaze_plot_tileset_proc_t *proc = &tmaze_ba_plot_pattern_tileset;
    tmaze_plot_pattern(txtcat("out/", dir), name, &M, TRUE, TRUE, proc, FALSE, 8.0, 1.0);
  }
  
void tmaze_figs_plt_ba_component_pattern(char *dir, char *prefix, int k, int which)
  {
    tmaze_pattern_t M;
    switch(k)
      {
      case 0: M = tmaze_figs_make_ba_component_pattern_0(which); break;
      case 1: M = tmaze_figs_make_ba_component_pattern_1(which); break;
      case 2: M = tmaze_figs_make_ba_component_pattern_2(which); break;
      case 3: M = tmaze_figs_make_ba_component_pattern_3(which); break;
      case 4: M = tmaze_figs_make_ba_component_pattern_4(which); break;
      default: 
        fprintf(stderr, "** %s not implemented for k = %d\n", __FUNCTION__, k);
        exit(1);
      }
      
    char *name = jsprintf("%s-%02d-%02d", prefix, k, which);
    tmaze_plot_tileset_proc_t *proc = &tmaze_ba_plot_pattern_tileset;
    tmaze_plot_pattern(txtcat("out/", dir), name, &M, TRUE, TRUE, proc, FALSE, 8.0, 1.0);
    free(name);
  }
  
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_0(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_ba_pattern_from_string("S;N;", FALSE, TRUE, 1); break;
      case 1: M = tmaze_ba_pattern_from_string("EW;",  FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }
  
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_1(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_ba_pattern_from_string("ENW;*N*;",  FALSE, TRUE, 1); break;
      case 1: M = tmaze_ba_pattern_from_string("S*;WW;N*;", FALSE, TRUE, 1); break;
      case 2: M = tmaze_ba_pattern_from_string("*S*;ESW;",  FALSE, TRUE, 1); break;
      case 3: M = tmaze_ba_pattern_from_string("*S;EE;*N;", FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }
  
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_2(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_ba_pattern_from_string("SS;WE;NN;",       FALSE, TRUE, 1); break;
      case 1: M = tmaze_ba_pattern_from_string("**S;ENE;*NN;",    FALSE, TRUE, 1); break;
      case 2: M = tmaze_ba_pattern_from_string("ENNW;*NN*;",      FALSE, TRUE, 1); break;
      case 3: M = tmaze_ba_pattern_from_string("**S*;ENSW;*N**;", FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }

tmaze_pattern_t tmaze_figs_make_ba_component_pattern_3(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 3: M = tmaze_ba_pattern_from_string("ENNW;*xE*;**N*;", FALSE, TRUE, 1); break;
      case 7: M = tmaze_ba_pattern_from_string("**S*;**WW;ENE*;*NN*;", FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d not implemented yet\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }

tmaze_pattern_t tmaze_figs_make_ba_component_pattern_4(int which)
  {
    tmaze_pattern_t M;
    switch(which)
      {
      case 0: M = tmaze_ba_pattern_from_string("ENNW;ESSW;",       FALSE, TRUE, 1); break;
      case 1: M = tmaze_ba_pattern_from_string("SS*;WE*;WWW;NN*;",       FALSE, TRUE, 1); break;
      default:
        fprintf(stderr, "** %s: pattern %d is undefined\n", __FUNCTION__, which);
        exit(1);
      }
    return M;
  }
