/* Drawing of asexual evolution results. */
#ifndef gdr_draw_H
#define gdr_draw_H
/* Last edited on 2023-06-25 19:03:41 by stolfi */

#define gdr_draw_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h> 

#include <gdr_demo.h>
#include <gdr_sim.h>

void gdr_draw_state_plot
  ( char *prefix,
    gdr_sim_state_t *st,
    int32_t y0,
    int32_t y1,
    int32_t yRef
  );
  /* Writes to file "{prefix}-evol.eps" an encapsulated
    Postscript figure depicting the result of {gdr_sim_run},
    clipped to the year range {y0..y1}. 
    
    The parameter {st} must be the state at the end of that simulation.
    If {yRef} is positive, draws the surviving lineages that count for 
    that year with distinctive colors. */
    
#define gdr_plot_state_file_INFO \
  "An Encapsulated Postscrio (EPS) file containing a plot of the evolution, showing" \
  " states as fat horizontal lines and parenting relations by arrows. The" \
  " horizontal axis is year, and each individual is shown as" \
  " a fat line from its birt year to its death year."

#endif

