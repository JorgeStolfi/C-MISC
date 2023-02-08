/* Last edited on 2012-12-08 23:37:21 by stolfilocal */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <zf.h>
#include <pswr.h>
#include <frgb.h>
#include <affirm.h>
/* #include <aa_trapez.h> */
/* #include <r3.h> */
/* #include <hr3.h> */
/* #include <ia.h> */
/* #include <ia_butfly.h> */
/* #include <ia_trapez.h> */
/* #include <aa.h> */
/* #include <flt.h> */
/* #include <pcode.h> */
/* #include <iaeval.h> */
/* #include <aaeval.h> */
/* #include <flteval.h> */
/* #include <jsstring.h> */
/* #include <jsfile.h> */

#include <irt_inter.h>

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

void irt_plot_interval(PSStream *ps, Interval *td, Interval *fd, Interval xv, Interval yv, zf_kind_t tv);
  /* Plots the rectangle {xv × yv} to {ps}, assuming that the overall plot area is 
     {td × fd}. The color is chosen according to the classification {tv}. */

/*** IMPLEMENTATIONS ***/

int irt_num_rays = 0;

void irt_compute_intersection
  ( shape_t *sh,         /* The object's shape */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *org,
    hr3_point_t *dst,
    Interval *hit,
    int *slo, int *shi,
    bool_t print_ray,
    PSStream *ps
  )
  {
    #define EPSILON 0.0
    #define DELTA 1.0e-6

    Interval unit_int = (Interval){Zero, One};
    zf_kind_t tn; /* Type of next interval after root */

    if (print_ray)
      { hr3_point_print(stderr, "  org = ", org, "%8.4f", "\n");
        hr3_point_print(stderr, "  dst = ", dst, "%8.4f", "\n");
      }

    irt_num_rays++;

    auto void eval_shape(Interval *tv, Interval *y, ia_butfly_t *a);
      /* Passed to {zf_enum_zeros} as the {eval} procedure.
        Evaluates the homogeneous four-variable function {sh} for all
        the points on the segment {data->org} to {data->dst} that correspond
        to {t} parameter values in the interval {x}.  

        This computation is performed with the arithmetic model
        {data->arith}.

        Returns in {y} the computed range of {f} in the interval {x}.
        Also stores in {a} a {ia_butfly_t} that contains the graph of
        {data->proc} along that ray segment. */

    auto bool_t process_interval (Interval *xv, Interval *yv, zf_kind_t tv);
      /*  Passed to {zf_enum_zeros} as the {report} procedure.  If the interval
        {xv} is a proper root, stores it in {data->hit} and tells {zf_enum_zeros}
        to stop looking for more roots.  Otherwise updates {data->slo}. */

    /* Overall ranges for plotting: */
    Interval td = (Interval) { Zero, One }; /* True range of {t} parameter. */
    Interval fd = (Interval) { +1, -1 };  /* Guessed range of {f}; to be set later. */

    /* initalize {hit} with empty interval: */
    hit->lo = One;
    hit->hi = Zero;

    /* Look for first root: */
    *slo = 0;
    tn = zf_enum_zeros(
      eval_shape,
      process_interval,
      unit_int,
      EPSILON,
      DELTA
    );

    /* Compute {shi}: */
    if (hit->lo <= hit->hi)
      { /* Found a proper root */
        if (tn == zf_kind_positive)
          { *shi = +1; }
        else if (tn == zf_kind_negative)
          { *shi = -1; }
        else if (tn == zf_kind_root)
          { fatalerror("irt_compute_intersection: contiguous roots"); }
        else
          { fatalerror("irt_compute_intersection: invalid next interval type"); }
      }
    else
      { /* No proper roots; */
        *shi = *slo;
      }

    if (print_ray)
      { fprintf(stderr, "  hit = (%+d) ", *slo);
        ia_print(stderr, *hit);
        fprintf(stderr, " (%+d)\n", *shi);
        fprintf(stderr, "\n");
      }
    return;

    #undef EPSILON
    #undef DELTA

    /* Internal functions: */
    
    bool_t process_interval(Interval *xv, Interval *yv, zf_kind_t tv)
      { if (ps != NULL) 
          { /* Plot the interval on {ps} */
            if (fd.lo > fd.hi)
              { /* Hack: set the function plot range {fd} from the first interval. */
                Float ff = FMAX(FABS(yv->lo), FABS(yv->hi));
                fd = (Interval){ -1.5*ff, +1.5*ff }; 
              }
            irt_plot_interval(ps, &td, &fd, *xv, *yv, tv);
          }
        affirm(hit->lo > hit->hi, "irt_process_interval: zf_enum_zeros didn't stop!");
        /* Set {sv} to +1 if interval is positive, -1 if negative, 0 if root-like: */
        int sv;
        switch (tv)
          { case zf_kind_positive: sv = +1; break;
            case zf_kind_negative: sv = -1; break;
            case zf_kind_root:     sv = 0; break;
            default: 
              fprintf(stderr, "irt_process_interval: invalid interval type");
              sv = 0;
          }
        if (print_ray)
          { fputc('\n', stderr);
            fprintf(stderr, "----------------------------------------\n");
            fprintf(stderr, "process_interval kind = %+2d:\n", sv);
            fputs("  t  =  ", stderr); ia_print(stderr, *xv); fputc('\n', stderr);
            fputs("  f  =  ", stderr); ia_print(stderr, *yv); fputc('\n', stderr);
            fprintf(stderr, "----------------------------------------\n");
            fputc('\n', stderr);
          }
        if (sv != 0)
          { *(slo) = sv; return(FALSE); }
        else if ((xv->lo > Zero) && (xv->hi < One))
          { *(hit) = *xv; return(TRUE); }
        else
          { return(FALSE); }
      }

    void eval_shape(Interval *tv, Interval *y, ia_butfly_t *a)
      { irt_num_evals++;
        irt_eval_shape_on_sub_segment(sh, arith, tv, print_ray, org, dst, y, a);
      }
  }

void irt_plot_interval(PSStream *ps, Interval *td, Interval *fd, Interval xv, Interval yv, zf_kind_t tv)
  {
    double xp[4], yp[4];

    /* Clip the rectangel to the Y range: */
    yv = ia_meet(yv, *fd);

    /* Plot the rectangle, scaling it to the {[0_1]} square: */
    ROUND_NEAR;
    double sx = 1.0/(td->hi - td->lo);
    double sy = 1.0/(fd->hi - fd->lo);
    xp[0] = sx*(xv.lo - td->lo);
    yp[0] = sy*(yv.lo - fd->lo);

    xp[1] = sx*(xv.hi - td->lo);
    yp[1] = sy*(yv.lo - fd->lo);

    xp[2] = sx*(xv.hi - td->lo);
    yp[2] = sy*(yv.hi - fd->lo);

    xp[3] = sx*(xv.lo - td->lo);
    yp[3] = sy*(yv.hi - fd->lo);

    /* Choose the color based on the interval's classification {tv}:  */
    frgb_t color;
    if (tv == zf_kind_positive)
      { color = (frgb_t){{ 1.00, 0.40, 0.20 }}; }
    else if (tv == zf_kind_negative)
      { color = (frgb_t){{ 0.20, 0.60, 1.00 }}; }
    else if (tv == zf_kind_undefined)
      { color = (frgb_t){{ 0.40, 0.40, 0.40 }}; }
    else if (tv == zf_kind_mixed)
      { color = (frgb_t){{ 0.50, 0.80, 0.10 }}; }
    else
      { color = (frgb_t){{ 1.00, 0.00, 0.20 }}; }

    pswr_set_fill_color(ps,  color.c[0], color.c[1], color.c[2]);
    pswr_polygon(ps, TRUE, xp, yp, 4, TRUE, TRUE, FALSE);
  }
