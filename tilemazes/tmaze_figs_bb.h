#ifndef tmaze_figs_bb_H
#define tmaze_figs_bb_H

/* Figures for the tile mazes paper - Blip-Blop mazes. */
/* Last edited on 2023-02-04 06:58:14 by stolfi */

#include <stdint.h>
#include <limits.h>
#include <stdio.h>

#include <bool.h>

#include <tmaze.h>
#include <tmaze_bb.h>

void tmaze_figs_plt_bb_tileset(char *dir, char *name, tmaze_tileset_t tset);
void tmaze_figs_plt_bb_sample(char *dir, char *name, int nx, int ny, int seed, bool_t showbig);
void tmaze_figs_plt_bb_conjuntive_pat(char *dir, char *name);
void tmaze_figs_plt_bb_component_pattern(char *dir, char *prefix, int k, int which);
  /* Each of these procedures writes an EPS file called "out/{dir}/{name}.eps"
    with a figure for the tile mazes paper (Blip Blop maze family).
    The {showbig} parameter highlights the large components. */
 
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_04(int which);
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_08(int which);
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_12(int which);
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_16(int which);
tmaze_pattern_t tmaze_figs_make_bb_component_pattern_20(int which);
  /* The procedure {tmaze_figs_make_bb_component_pattern_{k}} creates
    a patern that recognizes connected components of the Blip Blop
    maze of a specific T-shape, using {k} edges (but possibly covering
    less than {k} tiles). The t-shape is selected by the index
    {which}, starting from 0. Sets {M.torus} to FALSE, {M.sub} to
    TRUE, and {M.val} to 0. */
 
#endif
