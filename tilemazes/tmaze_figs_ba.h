#ifndef tmaze_figs_ba_H
#define tmaze_figs_ba_H

/* Figures for the tile mazes paper - Blip-Blop mazes. */
/* Last edited on 2023-02-04 06:57:55 by stolfi */

#include <stdint.h>
#include <limits.h>
#include <stdio.h>

#include <bool.h>

#include <tmaze.h>
#include <tmaze_ba.h>

void tmaze_figs_plt_ba_tileset(char *dir, char *name, tmaze_tileset_t tset);
void tmaze_figs_plt_ba_sample(char *dir, char *name, int nx, int ny, int seed, bool_t showbig);
void tmaze_figs_plt_ba_conjuntive_pat(char *dir, char *name);
void tmaze_figs_plt_ba_component_pattern(char *dir, char *prefix, int k, int which);
  /* Each of these procedures writes an EPS file called "out/{dir}/{name}.eps"
    with a figure for the tile mazes paper (Brasilia Airport maze family).
    The {showbig} parameter highlights the large components. */
 
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_0(int which);
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_1(int which);
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_2(int which);
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_3(int which);
tmaze_pattern_t tmaze_figs_make_ba_component_pattern_4(int which);
  /* The procedure {tmaze_figs_make_ba_component_pattern_{k}} creates
    a patern that recognizes connected components of the Brasilia
    Airport maze of a specific T-shape, covering {k} tiles. The
    t-shape is selected by the index {which}, starting from 0. Sets
    {M.torus} to FALSE, {M.sub} to TRUE, and {M.val} to 0. */
 
#endif
