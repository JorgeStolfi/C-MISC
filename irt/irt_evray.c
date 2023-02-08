/* Last edited on 2012-12-08 23:37:17 by stolfilocal */
/* See {irt_evray.h}. */

#include <irt_evray.h>

#include <zf.h>
#include <r3.h>
#include <hr3.h>
#include <ia.h>
#include <ia_butfly.h>
#include <ia_trapez.h>
#include <aa.h>
#include <aa_trapez.h>
#include <flt.h>
#include <pcode.h>
#include <iaeval.h>
#include <aaeval.h>
#include <flteval.h>
#include <pswr.h>
#include <affirm.h>
#include <jsstring.h>
/* #include <jsfile.h> */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

void irt_check_seg_eval_aa
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers.  */
    AAP *stack,           /* Evaluation stack  */
    hr3_point_t *org,     /* Start of ray */
    hr3_point_t *dst,     /* End of ray */
    AAP aat,              /* The affine form of the {t} parameter */
    ia_butfly_t *ip,      /* Bounding butterfly for f(aat) */
    AAP f                 /* The affine form for f(aat)  */
  );
  /* Checks the consistency of the affine form stack[0], computed inside
    {irt_seg_eval_aa}. */

/*** IMPLEMENTATIONS ***/

int irt_num_evals = 0;

Float irt_eval_f_on_point
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    Float *regs,          /* Evaluation registers */
    Float *stack,         /* Evaluation stack */
    hr3_point_t *org,     /* Origin of ray */
    hr3_point_t *dst,     /* End of ray */
    Float t               /* Parameter value where to evaluate */
  )
  { 
    double tt = (double) t;
    int i;

    ROUND_NEAR;

    hr3_point_t p = hr3_point_mix(1.0 - tt, org, tt, dst);
    for (i=0; i<4; i++) regs[i] = p.c.c[i];
    flt_eval(regs, stack, proc->code);
    return(stack[0]);
  }

void irt_eval_shape_on_sub_segment
  ( shape_t *sh, 
    arith_t arith, 
    Interval *tv, 
    bool_t print_ray,
    hr3_point_t *org,
    hr3_point_t *dst,
    Interval *y, 
    ia_butfly_t *a
  )
  { if (print_ray)
      { fprintf(stderr, "    t  =  "); ia_print(stderr, *tv);
        fprintf(stderr, "\n");
      }

    switch(arith)
      { case arith_ia:
          irt_seg_eval_ia (
            &(sh->proc),
            sh->ia_regs, sh->ia_stack,
            org, dst, tv, print_ray,
            a, y
          );
          break;
        case arith_ia_diff:
          irt_seg_eval_ia_diff (
            &(sh->proc),
            sh->id_regs, sh->id_stack,
            org, dst, tv, print_ray,
            a, y
          );
          break;
        case arith_aa:
          irt_seg_eval_aa (
            &(sh->proc),
            sh->aa_regs, sh->aa_stack,
            org, dst, tv, print_ray,
           a, y
          );
          break;
        case arith_mix:
          irt_seg_eval_mix (
            &(sh->proc),
            sh->ia_regs, sh->ia_stack,
            sh->aa_regs, sh->aa_stack,
            org, dst, tv, print_ray,
            a, y
          );
          break;
        default:
          fatalerror("irt_eval_shape_on_sub_segment: bad arithmetic");
      }

    if (print_ray) 
      { fprintf(stderr, "    f  =  ");
        ia_trapez_print(stderr, &(a->tp[0]));
        fprintf(stderr, "\n");
        fprintf(stderr, "    f  =  ");
        ia_trapez_print(stderr, &(a->tp[1]));
        fprintf(stderr, "\n");
        fprintf(stderr, "\n");
      }
  }

void irt_seg_eval_ia
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    Interval *regs,       /* Evaluation registers */
    Interval *stack,      /* Evaluation stack */
    hr3_point_t *org,     /* Start of ray */
    hr3_point_t *dst,     /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    ia_butfly_t *ip,      /* Out: bounding butterfly for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  )
  { int i;
    double a, b;
    for (i=0; i<4; i++)
      { a = org->c.c[i];
        b = dst->c.c[i];
        if (a < b)
          { ROUND_DOWN;
            regs[i].lo = a + tv->lo * (b - a);
            ROUND_UP;
            regs[i].hi = a + tv->hi * (b - a);
          }
        else
          { ROUND_DOWN;
            regs[i].lo = b + (1.0 - tv->hi) * (a - b);
            ROUND_UP;
            regs[i].hi = b + (1.0 - tv->lo) * (a - b);
          }
      }
    ia_eval (regs, stack, proc->code);
    (*ip) = ia_butfly_from_box(tv, &(stack[0])); 
    *yr = stack[0];
    ROUND_NEAR;
  }

void irt_seg_eval_ia_diff
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    IntervalDiff *regs,   /* Evaluation registers */
    IntervalDiff *stack,  /* Evaluation stack */
    hr3_point_t *org,      /* Start of ray */
    hr3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    ia_butfly_t *ip,      /* Out: bounding butterfly for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  )
  { int i;
    double a, b;
    /*   Note: the derivatives are relative to a new paramater {u}, which ranges from
      0 to 1 as {t} ranges from {tv->lo} to {tv->hi}. */
    for (i=0; i<4; i++)
      { a = org->c.c[i];
        b = dst->c.c[i];
        if (a < b)
          { ROUND_DOWN;
            regs[i].f.lo = a + tv->lo * (b - a);
            ROUND_UP;
            regs[i].f.hi = a + tv->hi * (b - a);
          }
        else
          { ROUND_DOWN;
            regs[i].f.lo = b + (1.0 - tv->hi) * (a - b);
            ROUND_UP;
            regs[i].f.hi = b + (1.0 - tv->lo) * (a - b);
          }
        ROUND_DOWN;
        regs[i].df.lo = (tv->hi - tv->lo)*(b - a);
        ROUND_UP;
        regs[i].df.hi = (tv->hi - tv->lo)*(b - a);
      }
    /* !!! Must make sure that the derivative is undef if {sh} is not differentiable !!! */
    ia_eval_diff (regs, stack, proc->code);
    *yr = stack[0].f;

    /* Compute bounding butterfly from function and derivative ranges: */
    (*ip) = ia_butfly_from_ia_diff(tv, ia_mid(*tv), &(stack[0].f), &(stack[0].df));
    ROUND_NEAR;
  }

void irt_seg_eval_aa
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack,           /* Evaluation stack */
    hr3_point_t *org,      /* Start of ray */
    hr3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    ia_butfly_t *ip,      /* Out: bounding butterfly for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  )
  { MemP frame = aa_top();
    AAP aat = aa_from_interval (*tv);
    AAP aar = aa_affine (aat, -1.0, 1.0, 1.0, 0.0); /* r = 1-t */
    AAP aaf;
    int i;

    if (print_ray)
      { fprintf(stderr, "      aat = "); aa_print(stderr, aat);
        fprintf(stderr, "\n");
        fprintf(stderr, "      aar = "); aa_print(stderr, aar);
        fprintf(stderr, "\n");
      }

    for (i=0; i<4; i++)
      { regs[i] = aa_affine_2 (aar, org->c.c[i], aat, dst->c.c[i], 1.0, 0.0, 0.0);
        if (print_ray)
          { fprintf(stderr, "      regs[%d] = ", i); aa_print(stderr, regs[i]);
            fprintf(stderr, "\n");
          }
      }
    aa_eval (regs, stack, proc->code);
    aaf = stack[0];
    *yr = aa_range(aaf);
    
    ia_trapez_t it = aa_trapez_from_pair(tv, aat, aaf);
    (*ip) = ia_butfly_from_trapez(&it);

    if (print_ray)
      { fprintf(stderr, "      result =    "); aa_print(stderr, stack[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "      butterfly = ");
        ia_butfly_print(stderr, ip, "  ");
        fprintf(stderr, "\n");
      }

    /*   if (print_ray && (aat->nterms != 0))
        { irt_check_seg_eval_aa(proc, regs, stack, org, dst, aat, ip, stack[0]); }
    */

    aa_flush(frame);
    ROUND_NEAR;
  }

void irt_check_seg_eval_aa
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack,           /* Evaluation stack */
    hr3_point_t *org,      /* Start of ray */
    hr3_point_t *dst,      /* End of ray */
    AAP aat,              /* The affine form of the {t} parameter */
    ia_butfly_t *ip,      /* Bounding butterfly for {f(aat)} */
    AAP f                 /* The affine form for {f(aat)}  */
  )
  {
    #define NCHECK 11

    MemP frame = aa_top();
    AATerm eps[1];
    int i, j;

    affirm((aat->nterms == 1), "irt_check_seg_eval_aa: bad aat");
    eps[0].id = ((AATermP)(aat + 1))->id;

    fprintf(stderr, "\n");
    for (j=-NCHECK; j<=NCHECK; j++)
      { /* The following comments assume that {t(e)} is the exact value of {t}
          corresponding to {eps_k = e}, where {eps_k} is the (unique)
          noise var that appears in {aat}. */

        /* Pick a sample value {e} in {[-1_+1]} for {eps_k}: */
        ROUND_NEAR; 
        Float e = ((Float) j)/((Float) NCHECK);
        
        /* Compute the affine form {ta} for the time {t(e)}: */
        eps[0].coef = e;
        AAP ta = aa_fix_eps(aat, 1, eps);
        
        /* Compute the affine form {sa} for {1-t(e)}: */
        AAP sa = aa_affine(ta, -1.0, 1.0, 1.0, 0.0);
        
        /* Compute the range {f_fix} for {f(t(e))} as {aa_fix_eps(f(aat), e)}: */
        Interval f_fix = aa_range(aa_fix_eps(f, 1, eps));
        
        /* Compute the range {f_cmp}for {f(t(e))} as {f(aa_fix_eps(aat, e))}: */ 
        for (i=0; i<4; i++)
          { regs[i] = aa_affine_2(sa, org->c.c[i], ta, dst->c.c[i], 1.0, 0.0, 0.0); }
        aa_eval (regs, stack, proc->code);
        Interval f_cmp = aa_range(stack[0]);
        
        /* Compute the range for {f(t(e))} obtained by slicing {ip} at {t(e)}: */
        Interval tev = aa_range(ta); /* A Small interval containing {t(e)}. */
        Interval f_trp = (Interval){ +1, -1 }; /* Enclosing interval for {f(t(e))}. */  
        int j;
        for (j = 0; j < 2; j++)
          { /* Get the trapezoid {j} from the butterfly, clipped to domain {tev}: */
            ia_trapez_t tpj = ia_trapez_clip(&tev, &(ip->tp[j]));
            if (tpj.x.lo <= tpj.x.hi)
              { /* The trapezoid {tpj} is not empty: */
                assert(tpj.x.lo == tev.lo);
                assert(tpj.x.hi == tev.hi);
                f_trp = ia_join(f_trp, tpj.yxlo);
                f_trp = ia_join(f_trp, tpj.yxhi);
              }
          }
            
        fprintf(stderr, "      ");
        fprintf(stderr, "  eps = "); ROUND_NEAR; flt_print(stderr, e);
        fprintf(stderr, "  t = "); ia_print(stderr, aa_range(ta));
        fprintf(stderr, "\n");
        fprintf(stderr, "      ");
        fprintf(stderr, "  f_cmp = "); ia_print(stderr, f_cmp);
        fprintf(stderr, "\n");
        fprintf(stderr, "      ");
        fprintf(stderr, "  f_fix = "); ia_print(stderr, f_fix);
        fprintf(stderr, "\n");
        fprintf(stderr, "      ");
        fprintf(stderr, "  f_trp = "); ia_print(stderr, f_trp);
        fprintf(stderr, "\n");
        fprintf(stderr, "\n");
        if ((f_cmp.lo > f_fix.hi) || (f_cmp.hi < f_fix.lo))
          { fatalerror("irt_check_seg_eval_aa: inconsistent f_cmp, f_fix"); }
        if ((f_trp.lo > f_fix.hi) || (f_trp.hi < f_fix.lo))
          { fatalerror("irt_check_seg_eval_aa: inconsistent f_trp, f_fix"); }
        if ((f_cmp.lo > f_trp.hi) || (f_cmp.hi < f_trp.lo))
          { fatalerror("irt_check_seg_eval_aa: inconsistent f_cmp, f_trp"); }
        aa_flush(frame);
      }
    fprintf(stderr, "\n");

    #undef NCHECK
  }

void irt_seg_eval_mix
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    Interval *iregs,      /* Evaluation registers (IA) */
    Interval *istack,     /* Evaluation stack (IA) */
    AAP *aregs,           /* Evaluation registers (AA) */
    AAP *astack,          /* Evaluation stack (AA) */
    hr3_point_t *org,      /* Start of ray */
    hr3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    ia_butfly_t *ip,      /* Out: bounding butterfly for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  )
  { /* Try first with IA, then use AA if undecided */
    irt_seg_eval_ia (proc, iregs, istack, org, dst, tv, print_ray, ip, yr);
    if( ((ip->tp[0].yxlo.lo <= Zero) || (ip->tp[1].yxhi.lo <= Zero)) &&
        ((ip->tp[0].yxlo.hi >= Zero) || (ip->tp[1].yxhi.hi >= Zero))
      )
      { /* The IA-computed range contains zero, recompute it with AA: */
        Interval yt;
        irt_seg_eval_aa (proc, aregs, astack, org, dst, tv, print_ray, ip, &yt);
        *yr = ia_meet(*yr, yt);
      }
  }

void irt_compute_surface_normal
  ( shape_t *sh,         /* The object's shape */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *hit,     /* The hit point. */
    r3_t *nrm            /* Out: The normal vector at {hit}. */
  )
  { int i;
    FloatDiff4 *regs  = sh->nrm_regs;
    FloatDiff4 *stack = sh->nrm_stack;

    for (i=0; i<4; i++)
      { FloatDiff4 *ri = &(regs[i]);
        ri->f = hit->c.c[i];
        ri->df[0] = 0.0;
        ri->df[1] = 0.0;
        ri->df[2] = 0.0;
        ri->df[3] = 0.0;
        ri->df[i] = 1.0;
      }

    flt_eval_diff4(regs, stack, sh->proc.code);

    nrm->c[0] = stack[0].df[1];
    nrm->c[1] = stack[0].df[2];
    nrm->c[2] = stack[0].df[3];
    (void) r3_dir(nrm, nrm);
    /* If {tgp} has the right orientation, {nrm} now points out. */

  }

void irt_debug_ray_graphically
  ( PSStream *ps,
    shape_t *sh,         /* The shape function */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *org,     /* Ray's origin */
    hr3_point_t *dst      /* Ray's destination */
  )
  {
    #define NUMINTS 32
    #define NUMSTEPS 256

    fprintf(stderr, "=== begin irt_debug_ray_graphically ===\n");
    
    /* Compute butterfiles {fvp[0..NUMINTS-1]} covering {F}, and its range {fd}: */
    ia_butfly_t fvp[NUMINTS];
    Interval tv[NUMINTS];

    Interval td = (Interval){Zero, One}; /* Parameter range. */
    Interval fd; /* Overall function range. */

    irt_eval_shape_on_equal_segments
      ( sh, arith, org, dst,
        NUMINTS, td, &fd, tv, fvp,
        TRUE
      );
  
    /* Widen the y plotting range a bit: */
    { Float r = Half * fd.hi - Half * fd.lo;
      fd.lo -= 0.05 * r;
      fd.hi += 0.05 * r;
    }
    if ((fd.lo <= MinusInfinity) || (fd.hi >= PlusInfinity))
      { fprintf(stderr, "irt_debug_ray_graphically: y range is infinite, aborted");
        return;
      }

    /* Start a new plot: */
    pswr_new_picture(ps, Zero, One, Zero, One);
    pswr_set_grid(ps, NUMINTS, 1);

    /* Now plot the butterfiles: */
    irt_plot_butterflies_along_segment(ps, NUMINTS, td, fd, tv, fvp);

    /* Draw the axes: */
    if ((td.lo < Zero) && (td.hi > Zero))
      { pswr_coord_line(ps, HOR, (Zero - td.lo)/(td.hi - td.lo)); }
    if ((fd.lo < Zero) && (fd.hi > Zero))
      { pswr_coord_line(ps, VER, (Zero - fd.lo)/(fd.hi - fd.lo)); }

    /* Now draw the function: */
    irt_plot_f_along_segment
      ( ps, sh, arith, org, dst,
        NUMSTEPS, td, fd
      );

    /* Draw the frame: */
    pswr_frame(ps);

    #undef NUMSTEPS
    #undef NUMINTS
    fprintf(stderr, "=== end irt_debug_ray_graphically ===\n");
  }

void irt_eval_shape_on_equal_segments
  ( shape_t *sh,         /* The shape function */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *org,    /* Ray origin */
    hr3_point_t *dst,    /* Ray destination */
    int nints,           /* Number of intervals to compute */
    Interval td,         /* Overall range of paramter (always [0_1] for now), */
    Interval *fd,        /* (Out) Overall range of function over {td} */
    Interval tv[],       /* (Out) Parameter intervals */
    ia_butfly_t fvp[],   /* (Out) Function butterfly for each parameter interval */
    bool_t print         /* TRUE prints the trapezoids. */
  )
  {
    int i;
    Interval yr;

    *fd = (Interval){Zero, Zero};

    for (i=0; i<nints; i++)
      { /* Compute the t interval: */
        if (i == nints-1)
          { tv[i].hi = td.hi; }
        else
          { ROUND_NEAR;
            tv[i].hi = td.lo + (td.hi - td.lo) * ((Float) (i+1))/((Float) nints);
          }

        if (i == 0)
          { tv[i].lo = td.lo; }
        else
          { tv[i].lo = tv[i-1].hi; }

        /* Evaluate F there: */
        irt_eval_shape_on_sub_segment
          ( sh, arith, &(tv[i]), FALSE, 
            org, dst, &yr, &(fvp[i])
          );
          
        if (print)
          { int j;
            fprintf(stderr, "  t[%02d]    = ", i);
            ia_print(stderr, tv[i]), 
            fprintf(stderr, "\n");
            for (j = 0; j < 2; j++)
              { fprintf(stderr, "  f[%02d][%d] =\n", i, j);
                ia_trapez_print(stderr, &(fvp[i].tp[j])); 
                fprintf(stderr, "\n");
              }
            fprintf(stderr, "\n");
          }

        /* Update the y range: */
        int j;
        for (j = 0; j < 2; j++)
          { if (! ia_is_full(&(fvp[i].tp[j].yxlo))) *fd = ia_join(*fd, fvp[i].tp[j].yxlo);
            if (! ia_is_full(&(fvp[i].tp[j].yxhi))) *fd = ia_join(*fd, fvp[i].tp[j].yxhi);
          }
      }

    if (fd->lo == fd->hi) *fd = (Interval){-Half, +Half};
  }

void irt_plot_butterflies_along_segment
  ( PSStream *ps,
    int nints,           /* Number of intervals to plot */
    Interval td,         /* Overall range of parameter */
    Interval fd,         /* Overall range of function over {td} */
    Interval tv[],       /* Parameter intervals (sub-intervals of {td}) */
    ia_butfly_t fvp[]    /* Function butterflies for each interval */
  )
  {
    double xp[4], yp[4];
    double gray = 0.75;
    int i;

    pswr_comment(ps, "begin plot of butterflies");

    for (i=0; i<nints; i++)
      {
        ia_butfly_t fvpi = fvp[i];
        Interval tvi = tv[i];
        assert(tvi.lo == fvpi.tp[0].x.lo);
        assert(tvi.hi == fvpi.tp[1].x.hi);
       
        int j;
        for (j = 0; j < 2; j++)
          { ia_trapez_t tpj = fvpi.tp[j]; 

            /* Scale trapezoid {tpj}  to [0 __ 1]: */
            if (ia_is_full(&(tpj.yxlo)) || ia_is_full(&(tpj.yxhi)))
              { /* If either end is infinite, replace the trapezoid by the clip box: */
                tpj.yxlo = fd; tpj.yxhi = fd;
              }

            /* Plot the trapezoid: */
            ROUND_NEAR;
            double sx = 1.0/(td.hi - td.lo);
            double sy = 1.0/(fd.hi - fd.lo);
            
            xp[0] = sx*(tpj.x.lo - td.lo);
            yp[0] = sy*(tpj.yxlo.lo - fd.lo);
            
            xp[1] = sx*(tpj.x.hi - td.lo);
            yp[1] = sy*(tpj.yxhi.lo - fd.lo);
            
            xp[2] = sx*(tpj.x.hi - td.lo);
            yp[2] = sy*(tpj.yxhi.hi - fd.lo);
            
            xp[3] = sx*(tpj.x.lo - td.lo);
            yp[3] = sy*(tpj.yxlo.hi - fd.lo);

            pswr_set_fill_color(ps,  gray,gray,gray);
            pswr_polygon(ps, TRUE, xp, yp, 4, TRUE, TRUE, FALSE);
          }
      }

    pswr_comment(ps, "end plot of butterflies");
  }

void irt_plot_f_along_segment
  ( PSStream *ps,
    shape_t *sh,         /* The shape function */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *org,     /* Ray origin */
    hr3_point_t *dst,     /* Ray destination */
    int nsteps,          /* Number of intervals to plot */
    Interval td,         /* Overall range of parameter */
    Interval fd          /* Overall range of function over {td} */
  )
  {
    Float t0, f0, t1, f1;
    double x0, y0, x1, y1;
    int i;

    ROUND_NEAR;
    pswr_comment(ps, "begin plot of actual graph");

    pswr_set_pen(ps, 0.00, 0.30, 0.00, 0.0,0.0,0.0);

    t1 = td.lo;
    f1 = irt_eval_f_on_point(&(sh->proc), sh->flt_regs, sh->flt_stack, org, dst, t1);

    for (i=0; i<nsteps; i++)
      {
        t0 = t1;
        f0 = f1;

        t1 = td.lo + ((td.hi - td.lo)*(i+1))/((Float) nsteps);
        f1 = irt_eval_f_on_point(&(sh->proc), sh->flt_regs, sh->flt_stack, org, dst, t1);

        if ( ((abs(f1) < Infinity) && (abs(f0) < Infinity))
          && ((f0 <= fd.hi) || (f1 <= fd.hi))
          && ((f0 >= fd.lo) || (f1 >= fd.lo))
          )
          { x0 = (t0 - td.lo)/(td.hi - td.lo); y0 = (f0 - fd.lo)/(fd.hi - fd.lo);
            x1 = (t1 - td.lo)/(td.hi - td.lo); y1 = (f1 - fd.lo)/(fd.hi - fd.lo);
            pswr_segment(ps, x0, y0, x1, y1);
          }

      }

    pswr_comment(ps, "end plot of actual graph");
  }

