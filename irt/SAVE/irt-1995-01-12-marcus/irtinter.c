#include "irtinter.h"
#include <zf.h>
#include <r3.h>
#include <h3.h>
#include <ia.h>
#include <aa.h>
#include <flt.h>
#include <pcode.h>
#include <iaeval.h>
#include <aaeval.h>
#include <flteval.h>
#include <ps.h>
#include <affirm.h>
#include <ioprotos.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

typedef struct {
    /* Parameters passed to irt_compute_intersection: */
    shape_t *sh;
    h3_point_t *org;
    h3_point_t *dst;
    Interval *hit;
    int *slo; int *shi;
    int print_ray;
  } cmp_int_data;

Trapezoid irt_eval_f_on_sub_segment (Interval *xv, void *data, Interval *yr);
  /*
    Passed to $zf_enum_zeros$ as the $f$ procedure.
    Evaluates the homogeneous four-variable function $data->sh$ for all
    the points on the segment $data->org$ to $data->dst$ that correspond
    to $t$ parameter values in the interval $xv$.

    Does so by computing the points $plo$ and $phi$ of the segment
    that correspond to $t=xv.lo$ and $t=xv.hi$, respectively.
    Then converts these two points into four coordinate intervals
    or affine forms, and calls the appropriate $proc$ evaluator.

    Returns an Trapezoid that describes a trapezoid that contains
    the graph of $data->proc$ along that ray segment and returns a rectangle $yr$
    that contains the range of $f$ at the interval $tv$. */

int irt_process_root (Interval *xv, Interval *yv, zf_type tv, void *data);
  /*
    Passed to $zf_enum_zeros$ as the $report$ procedure.  If the interval
    $xv$ is a proper root, stores it in $data->hit$ and tells zf_enum_zeros
    to stop looking for more roots.  Otherwise updates $data->slo$. */

void irt_trapezoid_from_ia_diff (IntervalDiff *fd, Trapezoid *ip);
  /*
    Computes a bounding trapezoid $ip$ from the ranges $fd->f$ and $fd->df$
    of a function and its derivativa.
    */

void irt_trapezoid_from_aa (Interval *tv, AAP aat, AAP aaf, Trapezoid *ip);
  /*
    Computes a bounding trapezoid $ip$ for a one-parameter function $f$,
    given a parameter interval $tv$, an affine form $aat$ that covers
    that parameter, and an affine form $aaf$ for $f(aat)$. */

Float irt_interpolate_up (
    Float xa, Float ya,  /* Leftmost point */
    Float xb, Float yb,  /* Rightmost point */
    Float x              /* Interpolation abscissa */
  );
  /*
    Interpolates between the points (xa, ya) and (xb,yb) at the abscissa x,
    rounding up. */

Float irt_interpolate_down (
    Float xa, Float ya,  /* Leftmost point */
    Float xb, Float yb,  /* Rightmost point */
    Float x              /* Interpolation abscissa */
  );
  /*
    Interpolates between the points (xa, ya) and (xb,yb) at the abscissa x,
    rounding down. */

void irt_seg_eval_ia (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    Interval *regs,       /* Evaluation registers */
    Interval *stack,      /* Evaluation stack */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    Trapezoid *ip,        /* Out: bounding trapezoid for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  );

void irt_seg_eval_ia_diff (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    IntervalDiff *regs,   /* Evaluation registers */
    IntervalDiff *stack,  /* Evaluation stack */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    Trapezoid *ip,        /* Out: bounding trapezoid for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  );

void irt_seg_eval_aa (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack,           /* Evaluation stack */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    Trapezoid *ip,        /* Out: bounding trapezoid for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  );

void irt_check_seg_eval_aa(
    pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers.  */
    AAP *stack,           /* Evaluation stack  */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    AAP aat,              /* The affine form of the $t$ parameter */
    Trapezoid *ip,        /* Bounding trapezoid for f(aat) */
    AAP f                 /* The affine form for f(aat)  */
  );
  /*
    Checks the consistency of the affine form stack[0], computed inside
    irt_seg_eval_aa. */

void irt_eval_f_on_equal_segments(
    shape_t *sh,         /* The shape function */
    h3_point_t *org,     /* Ray origin */
    h3_point_t *dst,     /* Ray destination */
    int nints,           /* Number of intervals to compute */
    Interval td,         /* Overall range of paramter (always [0_1] for now), */
    Interval *fd,        /* (Out) Overall range of function over $td$ */
    Interval tv[],       /* (Out) Parameter intervals */
    Trapezoid fvp[]      /* (Out) Function trapezoid for each parameter interval */
  );
  /*
    Divides the interval $td$ into $nints$ equal intervals, and
    evaluates $sh$ on the corresponding segments of the $org--dst$ ray,
    using $sh->arithmetic$. Returns the parameter segments in $tv[0..nints-1]$,
    the function trapezoids in $fvp[0..nints-1]$, and the overall
    range of the function in $fd$. */

void irt_plot_trapezoids_along_segment(
    FILE *psfile,
    int nints,           /* Number of intervals to plot */
    Interval td,         /* Overall range of parameter */
    Interval fd,         /* Overall range of function over $td$ */
    Interval tv[],       /* Parameter intervals (sub-intervals of $td$) */
    Trapezoid fvp[]      /* Function trapezoids for each interval */
  );
  /*
    Plots the trapezoids computed by irt_eval_f_on_equal_segments
    to the Postscript file $psfile$ */

void irt_plot_f_along_segment(
    FILE *psfile,
    shape_t *sh,         /* The shape function */
    h3_point_t *org,     /* Ray origin */
    h3_point_t *dst,     /* Ray destination */
    int nsteps,          /* Number of points to evaluate */
    Interval td,         /* Overall range of parameter */
    Interval fd          /* Overall range of function over $td$ */
  );
  /*
    Plots a graph of $sh->proc$ along the ray $org--dst$, evaluating it
    with ordinary floating-point on $nsteps$ points. */

Float irt_eval_f_on_point(
    pcode_proc_t *proc,   /* The function's pseudo-code */
    Float *regs,          /* Evaluation registers */
    Float *stack,         /* Evaluation stack */
    h3_point_t *org,      /* Origin of ray */
    h3_point_t *dst,      /* End of ray */
    Float t               /* Parameter value where to evaluate */
  );
  /*
    Evaluates $sh->proc$ at the point of ray $org--dst$ with parameter value $t$,
    using ordinary floating-point. */

/*** IMPLEMENTATIONS ***/

int irt_num_rays = 0;
int irt_num_evals = 0;

void irt_compute_intersection (
    shape_t *sh,         /* The object's shape */
    h3_point_t *org,
    h3_point_t *dst,
    Interval *hit,
    int *slo, int *shi,
    int print_ray
  )
  {
    #define EPSILON 0.0
    #define DELTA 1.0e-6

    cmp_int_data cd;
    Interval unit_int = (Interval){Zero, One};
    zf_type tn; /* Type of next interval after root */

    if (print_ray)
      { fprintf(stderr, "  org = "); h3_print_point(stderr, org);
        fprintf(stderr, "\n");
        fprintf(stderr, "  dst = "); h3_print_point(stderr, dst);
        fprintf(stderr, "\n");
      }

    irt_num_rays++;

    cd.sh = sh;
    cd.org = org;
    cd.dst = dst;
    cd.hit = hit;
    cd.slo = slo;
    cd.shi = shi;
    cd.print_ray = print_ray;

    /* initalize $hit$ with empty interval: */
    hit->lo = One;
    hit->hi = Zero;

    /* Look for first root: */
    *slo = 0;
    tn = zf_enum_zeros(
      irt_eval_f_on_sub_segment,
      irt_process_root,
      (void *) &cd,
      unit_int,
      EPSILON,
      DELTA
    );

    /* Compute $shi$: */
    if (hit->lo <= hit->hi)
      { /* Found a proper root */
        if (tn == zf_type_positive)
          { *shi = +1; }
        else if (tn == zf_type_negative)
          { *shi = -1; }
        else if (tn == zf_type_root)
          { error("irt_compute_intersection: contiguous roots"); }
        else
          { error("irt_compute_intersection: invalid next interval type"); }
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

    #undef EPSILON
    #undef DELTA

  }

int irt_process_root (Interval *xv, Interval *yv, zf_type tv, void *data)
  {
    cmp_int_data *cd = (cmp_int_data *) data;
    int sv;

    affirm(cd->hit->lo > cd->hit->hi, "irt_process_root: zf_enum_zeros didn't stop!");

    switch (tv)
      { case zf_type_positive: sv = +1; break;
        case zf_type_negative: sv = -1; break;
        case zf_type_root:     sv = 0; break;
        default: error("irt_process_root: invalid interval type");
      }

    if (sv != 0)
      { *(cd->slo) = sv; return(0); }
    else if ((xv->lo > Zero) && (xv->hi < One))
      { *(cd->hit) = *xv; return(1); }
    else
      { return(0); }
  }

Trapezoid irt_eval_f_on_sub_segment (Interval *xv, void *data, Interval *yr)
  { cmp_int_data *cd = (cmp_int_data *) data;
    Trapezoid ip;

    irt_num_evals++;

    if (cd->print_ray)
      { fprintf(stderr, "    xv =  "); ia_print(stderr, *xv);
        fprintf(stderr, "\n");
      }

    switch(cd->sh->arithmetic)
      {
        case arith_ia:
          irt_seg_eval_ia (
            &(cd->sh->proc), cd->sh->ia_regs, cd->sh->ia_stack,
            cd->org, cd->dst, xv, cd->print_ray,
            &ip, yr
          );
          break;
        case arith_ia_diff:
          irt_seg_eval_ia_diff (
            &(cd->sh->proc), cd->sh->id_regs, cd->sh->id_stack,
            cd->org, cd->dst, xv, cd->print_ray,
            &ip, yr
          );
          break;
        case arith_aa:
          irt_seg_eval_aa (
            &(cd->sh->proc), cd->sh->aa_regs, cd->sh->aa_stack,
            cd->org, cd->dst, xv, cd->print_ray,
           &ip, yr
          );
          break;
        default:
          error("irt_eval_f_on_sub_segment: bad arithmetic");
      }

    if (cd->print_ray)
      { fprintf(stderr, "      lox = ");
        ia_print(stderr, ip.lox);
        fprintf(stderr, "  hix = ");
        ia_print(stderr, ip.hix);
        fprintf(stderr, "\n");
      }

    return(ip);
  }

void irt_seg_eval_ia (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    Interval *regs,       /* Evaluation registers */
    Interval *stack,      /* Evaluation stack */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    Trapezoid *ip,        /* Out: bounding trapezoid for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  )
  { int i;
    double a, b;
    for (i=0; i<4; i++)
      { a = org->c[i];
        b = dst->c[i];
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
    ip->lox = ip->hix = stack[0];
    *yr = stack[0];
    ROUND_NEAR;
  }

void irt_seg_eval_ia_diff (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    IntervalDiff *regs,   /* Evaluation registers */
    IntervalDiff *stack,  /* Evaluation stack */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    Trapezoid *ip,        /* Out: bounding trapezoid for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  )
  { int i;
    double a, b;
    /*
      Note: the derivatives are relative to a new paramater $u$, which ranges from
      0 to 1 as $t$ ranges from $tv->lo$ to $tv->hi$. */
    for (i=0; i<4; i++)
      { a = org->c[i];
        b = dst->c[i];
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
    ia_eval_diff (regs, stack, proc->code);
    *yr = stack[0].f;

    /* Compute bounding trapezoid from function and derivative ranges: */
    irt_trapezoid_from_ia_diff(&(stack[0]), ip);
    ROUND_NEAR;
  }

void irt_trapezoid_from_ia_diff (IntervalDiff *fd, Trapezoid *ip)
  { ip->lox = ip->hix = fd->f;

    if (fd->df.lo > 0.0)
      { /* Function is strictly increasing */
        ROUND_UP;
        ip->lox.hi = ip->hix.hi - fd->df.lo;
        ROUND_DOWN;
        ip->hix.lo = ip->lox.lo + fd->df.lo;
      }
    else if (fd->df.hi < 0.0)
      { /* Function is strictly decreasing */
        ROUND_UP;
        ip->hix.hi = ip->lox.hi + fd->df.hi;
        ROUND_DOWN;
        ip->lox.lo = ip->hix.lo - fd->df.hi;
      }
    else
      { /* Function may be increasing or decreasing, can't say anything */ }
  }

void irt_seg_eval_aa (
    pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack,           /* Evaluation stack */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    Trapezoid *ip,        /* Out: bounding trapezoid for function along segment */
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
      { regs[i] = aa_affine_2 (aar, org->c[i], aat, dst->c[i], 1.0, 0.0, 0.0);
        if (print_ray)
          { fprintf(stderr, "      regs[%d] = ", i); aa_print(stderr, regs[i]);
            fprintf(stderr, "\n");
          }
      }
    aa_eval (regs, stack, proc->code);
    aaf = stack[0];
    *yr = aa_range(aaf);
    irt_trapezoid_from_aa (tv, aat, aaf, ip);

    if (print_ray)
      { fprintf(stderr, "      result =   "); aa_print(stderr, stack[0]);
        fprintf(stderr, "\n");
        fprintf(stderr, "      trapezoid = ");
        ia_print(stderr, ip->lox); fprintf(stderr, "  "); ia_print(stderr, ip->hix);
        fprintf(stderr, "\n");
      }

    /*
      if (print_ray && (aat->nterms != 0))
        { irt_check_seg_eval_aa(proc, regs, stack, org, dst, aat, ip, stack[0]); }
    */

    aa_flush(frame);
    ROUND_NEAR;
  }

void irt_check_seg_eval_aa(
    pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack,           /* Evaluation stack */
    h3_point_t *org,      /* Start of ray */
    h3_point_t *dst,      /* End of ray */
    AAP aat,              /* The affine form of the $t$ parameter */
    Trapezoid *ip,        /* Bounding trapezoid for f(aat) */
    AAP f                 /* The affine form for f(aat)  */
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
      { /*
          The following assumes that $t(e)$ is the exact value of $t$
          corresponding to $eps_k = e$, where $eps_k$ is the (unique)
          noise var that appears in $aat$.
        */
        Float e; /* A sample value for $eps_k$ */
        Interval f_cmp; /* Range for $f(t(e))$ computed as f(aa_fix_eps(aat, e)) */
        Interval f_fix; /* Range for $f(t(e))$ computed as aa_fix_eps(f(aat), e) */
        Interval f_trp; /* Range for $f(t(e))$ obtained by interpolat ing $ip$ */
        AAP ct; /* An affine form for $t(e)$ */
        AAP cr; /* An affine form for $1 - t(e)$. */

        ROUND_NEAR; e = ((Float) j)/((Float) NCHECK);
        eps[0].coef = e;

        ct = aa_fix_eps(aat, 1, eps);
        cr = aa_affine (ct, -1.0, 1.0, 1.0, 0.0);
        f_fix = aa_range(aa_fix_eps(f, 1, eps));
        for (i=0; i<4; i++)
          { regs[i] = aa_affine_2 (cr, org->c[i], ct, dst->c[i], 1.0, 0.0, 0.0); }
        aa_eval (regs, stack, proc->code);
        f_cmp = aa_range(stack[0]);
        { Float alo, ahi, blo, bhi;
          if (ip->lox.lo >= Zero)
            { ROUND_DOWN; alo = ((One - e)/Two)*ip->lox.lo;}
          else
            { ROUND_UP;   alo = -(((One - e)/Two)*(-ip->lox.lo)); }
          if (ip->hix.lo >= Zero)
            { ROUND_DOWN; blo = ((One + e)/Two)*ip->hix.lo;}
          else
            { ROUND_UP;   blo = -(((One + e)/Two)*(-ip->hix.lo)); }

          if (ip->lox.hi >= Zero)
            { ROUND_UP;   ahi = ((One - e)/Two)*ip->lox.hi;}
          else
            { ROUND_DOWN; ahi = -(((One - e)/Two)*(-ip->lox.hi)); }
          if (ip->hix.hi >= Zero)
            { ROUND_UP;   bhi = ((One + e)/Two)*ip->hix.hi;}
          else
            { ROUND_DOWN; bhi = -(((One + e)/Two)*(-ip->hix.hi)); }

          ROUND_DOWN; f_trp.lo = alo + blo;
          ROUND_UP;   f_trp.hi = ahi + bhi;
        }

        fprintf(stderr, "      ");
        fprintf(stderr, "  eps = "); ROUND_NEAR; flt_print(stderr, e);
        fprintf(stderr, "  t = "); ia_print(stderr, aa_range(ct));
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
          { error("irt_check_seg_eval_aa: inconsistent f_cmp, f_fix"); }
        if ((f_trp.lo > f_fix.hi) || (f_trp.hi < f_fix.lo))
          { error("irt_check_seg_eval_aa: inconsistent f_trp, f_fix"); }
        if ((f_cmp.lo > f_trp.hi) || (f_cmp.hi < f_trp.lo))
          { error("irt_check_seg_eval_aa: inconsistent f_cmp, f_trp"); }
        aa_flush(frame);
      }
    fprintf(stderr, "\n");

    #undef NCHECK
  }

void irt_trapezoid_from_aa (Interval *tv, AAP aat, AAP aaf, Trapezoid *ip)
  {
    if (aat->nterms == 0)
      { ip->lox = aa_implicit_range(aaf);
        ip->hix = ip->lox;
      }
    else
      {
        AATerm eps[1];
        Interval atv = aa_implicit_range(aat);
        /*
          Let RANGE(aa) be the range of an affine form $aa$,
          and FIX(aa, k, v) the result of substituting $v$ for $e_k$ in aa,
          both computed with exact (infinite-precison) arithmetic.

          Let R' be the joint range of the affine forms $aat$ and
          $aaf$.  Since $aat$ depends on only one epsilon $e_i$, this
          region is actually a trapezoid with vertical bases which
          contains the graph of $f(t)$ for $t$ in RANGE(aat). The
          projection of R' on the $t$-axis is RANGE(aat), and the two
          bases of R' are RANGE(FIX(aaf, i, -1)) and RANGE(FIX(aaf, i, +1)).

          We first compute a trapezoid R'' that contains R',
          and has the same projection on the $t$-axis.
          The bases of this trapezoid are computed by the formulas
          above, except that we use the conservative estimators
          aa_range and aa_fix_eps instead of RANGE and FIX.
        */
        eps[0].id = ((AATermP)(aat + 1))->id;
        eps[0].coef = -One;
        ip->lox = aa_implicit_range(aa_fix_eps(aaf, 1, eps));
        eps[0].coef = One;
        ip->hix = aa_implicit_range(aa_fix_eps(aaf, 1, eps));
        /*
          However, "zf_enum_zeros" expects a trapezoid whose
          $t$-projection is $tv$, not RANGE(aat).
          We must adjust $ip.lox$ and $ip.hix$ to match
          this expectation.
        */
        if (ip->lox.hi > ip->hix.hi)
          { ip->hix.hi = irt_interpolate_up(tv->lo, ip->lox.hi, atv.hi, ip->hix.hi, tv->hi); }
        else if (ip->lox.hi < ip->hix.hi)
          { ip->lox.hi = irt_interpolate_up(atv.lo, ip->lox.hi, tv->hi, ip->hix.hi, tv->lo); }
        else
          { /* ip->lox.hi, ip->hix.hi are  OK */ }
        if (ip->lox.lo < ip->hix.lo)
          { ip->hix.lo = irt_interpolate_down(tv->lo, ip->lox.lo, atv.hi, ip->hix.lo, tv->hi); }
        else if (ip->lox.lo > ip->hix.lo)
          { ip->lox.lo = irt_interpolate_down(atv.lo, ip->lox.lo, tv->hi, ip->hix.lo, tv->lo); }
        else
          { /* ip->lox.lo, ip->hix.lo are OK */ }
      }
  }

Float irt_interpolate_down (
    Float xa, Float ya,  /* Leftmost point */
    Float xb, Float yb,  /* Rightmost point */
    Float x              /* Interpolation abscissa */
  )
  {
    Float dx, dy, rx;
    affirm (xb > xa, "irt_interpolate_down: xb < xa!");
    ROUND_UP;
    dx = xb - xa;
    if (ya <= yb)
      { ROUND_DOWN;
        dy = yb - ya;
        rx = x - xa;
        return (ya + dy * (rx/dx));
      }
    else
      { ROUND_DOWN;
        dy = ya - yb;
        rx = xb - x;
        return (yb + dy * (rx/dx));
      }
 }

Float irt_interpolate_up (
    Float xa, Float ya,  /* Leftmost point */
    Float xb, Float yb,  /* Rightmost point */
    Float x              /* Interpolation abscissa */
  )
  {
    Float dx, dy, rx;
    affirm (xb > xa, "irt_interpolate_up: xb < xa!");
    ROUND_DOWN;
    dx = xb - xa;
    if (ya <= yb)
      { if (dx == Zero) return (yb);  /* Just in case... */
        ROUND_UP;
        dy = yb - ya;
        rx = x - xa;
        return (ya + dy * (rx/dx));
      }
    else
      { if (dx == Zero) return (ya);  /* Just in case... */
        ROUND_UP;
        dy = ya - yb;
        rx = xb - x;
        return (yb + dy * (rx/dx));
      }
  }

void irt_compute_surface_normal (
    shape_t *sh,         /* The object's shape */
    h3_point_t *hit,    /* The hit point. */
    r3_t *nrm           /* Out: The normal vector at $hit$. */
  )
  { int i;
    FloatDiff4 *regs  = sh->nrm_regs;
    FloatDiff4 *stack = sh->nrm_stack;

    for (i=0; i<4; i++)
      { FloatDiff4 *ri = &(regs[i]);
        ri->f = hit->c[i];
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
    (void) r3_normalize(nrm);
    /* If $tgp$ has the right orientation, $nrm$ now points out. */

  }

void irt_debug_ray_graphically (
    shape_t *sh,         /* The shape function */
    h3_point_t *org,     /* Ray's origin */
    h3_point_t *dst,     /* Ray's destination */
    char *scene_name,
    char *ray_name
  )
  {
    #define NUMINTS 32
    #define NUMSTEPS 256

    Trapezoid fvp[NUMINTS];
    Interval tv[NUMINTS];

    Interval td = (Interval){Zero, One};
    Interval fd;

    irt_eval_f_on_equal_segments(
      sh, org, dst,
      NUMINTS, td, &fd, tv, fvp
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

    /* Now generate the plot file: */

    {
      char *file_name = txtcat(txtcat(scene_name, "-"), txtcat(ray_name, ".ps"));

      FILE *psfile;

      double scale = 6.0 * 72.0;
      double xc = 4.25 * 72.0;
      double hmin = xc - scale / 2.0;
      double hmax = xc + scale / 2.0;
      double yc = 6.50 *72.0;
      double vmin = yc - scale / 2.0;
      double vmax = yc + scale / 2.0;

      psfile = open_write(file_name, TRUE);
      ps_begin_document(psfile, "letter");

      ps_begin_page(psfile, 1);
      ps_set_window(psfile,
        Zero, One,    Zero, One,
        hmin,  hmax,  vmin,  vmax,
        NUMINTS, 1
      );

      irt_plot_trapezoids_along_segment(psfile, NUMINTS, td, fd, tv, fvp);

      /* Draw the axes: */
      if ((td.lo < Zero) && (td.hi > Zero))
        { ps_draw_coord_line(psfile, 'x', (Zero - td.lo)/(td.hi - td.lo)); }
      if ((fd.lo < Zero) && (fd.hi > Zero))
        { ps_draw_coord_line(psfile, 'y', (Zero - fd.lo)/(fd.hi - fd.lo)); }

      /* Now draw the function: */
      irt_plot_f_along_segment (
        psfile, sh, org, dst,
        NUMSTEPS, td, fd
      );

      ps_draw_frame(psfile);
      ps_end_page(psfile);
      ps_end_document(psfile, 1);
      fclose(psfile);

    }

    #undef NUMSTEPS
    #undef NUMINTS

  }

void irt_eval_f_on_equal_segments(
    shape_t *sh,         /* The shape function */
    h3_point_t *org,     /* Ray origin */
    h3_point_t *dst,     /* Ray destination */
    int nints,           /* Number of intervals to compute */
    Interval td,         /* Overall range of paramter (always [0_1] for now), */
    Interval *fd,        /* (Out) Overall range of function over $td$ */
    Interval tv[],       /* (Out) Parameter intervals */
    Trapezoid fvp[]      /* (Out) Function trapezoid for each parameter interval */
  )
  {
    cmp_int_data cd;
    int i;
    Interval yr;

    cd.sh = sh;
    cd.org = org;
    cd.dst = dst;
    cd.hit = NULL;
    cd.slo = NULL;
    cd.shi = NULL;
    cd.print_ray = 0;

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
        fvp[i] = irt_eval_f_on_sub_segment(&(tv[i]), &cd, &yr);

        /* Update the y range: */
        if (! ia_is_full(&(fvp[i].lox))) *fd = ia_join(*fd, fvp[i].lox);
        if (! ia_is_full(&(fvp[i].hix))) *fd = ia_join(*fd, fvp[i].hix);
      }

    if (fd->lo == fd->hi) *fd = (Interval){-Half, +Half};
  }

void irt_plot_trapezoids_along_segment(
    FILE *psfile,
    int nints,           /* Number of intervals to plot */
    Interval td,         /* Overall range of parameter */
    Interval fd,         /* Overall range of function over $td$ */
    Interval tv[],       /* Parameter intervals (sub-intervals of $td$) */
    Trapezoid fvp[]      /* Function trapezoids for each interval */
  )
  {
    double xp[4], yp[4];
    double gray = 0.75;
    int i;

    for (i=0; i<nints; i++)
      {
        Trapezoid fvpi = fvp[i];
        Interval tvi = tv[i];

        /* Scale f trapezoid to [0 __ 1]: */
        if (ia_is_full(&fvpi.lox) || ia_is_full(&fvpi.hix))
          { fvpi.lox = fd; fvpi.hix = fd; }

        /* Plot trapezoid: */
        ROUND_NEAR;
        xp[0] = (tvi.lo - td.lo)/(td.hi - td.lo);
        yp[0] = (fvpi.lox.lo - fd.lo)/(fd.hi - fd.lo);
        xp[1] = (tvi.hi - td.lo)/(td.hi - td.lo);
        yp[1] = (fvpi.hix.lo - fd.lo)/(fd.hi - fd.lo);
        xp[2] = (tvi.hi - td.lo)/(td.hi - td.lo);
        yp[2] = (fvpi.hix.hi - fd.lo)/(fd.hi - fd.lo);
        xp[3] = (tvi.lo - td.lo)/(td.hi - td.lo);
        yp[3] = (fvpi.lox.hi - fd.lo)/(fd.hi - fd.lo);

        ps_fill_and_draw_polygon(psfile, xp, yp, 4, gray);
      }

    ps_end_section(psfile);
  }

void irt_plot_f_along_segment(
    FILE *psfile,
    shape_t *sh,         /* The shape function */
    h3_point_t *org,     /* Ray origin */
    h3_point_t *dst,     /* Ray destination */
    int nsteps,          /* Number of intervals to plot */
    Interval td,         /* Overall range of parameter */
    Interval fd          /* Overall range of function over $td$ */
  )
  {
    Float t0, f0, t1, f1;
    double x0, y0, x1, y1;
    int i;

    ROUND_NEAR;
    ps_begin_section(psfile, "Plot of actual graph");

    ps_set_pen(psfile, 0.0, 0.30, 0.0, 0.0);

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
            ps_draw_segment(psfile, x0, y0, x1, y1);
          }

      }

    ps_end_section(psfile);
  }

Float irt_eval_f_on_point(
    pcode_proc_t *proc,   /* The function's pseudo-code */
    Float *regs,          /* Evaluation registers */
    Float *stack,         /* Evaluation stack */
    h3_point_t *org,      /* Origin of ray */
    h3_point_t *dst,      /* End of ray */
    Float t               /* Parameter value where to evaluate */
  )
  { h3_point_t p;
    double tt = (double) t;
    int i;

    ROUND_NEAR;
    h3_mix (1.0 - tt, org, tt, dst, &p);
    for (i=0; i<4; i++) regs[i] = p.c[i];
    flt_eval(regs, stack, proc->code);
    return(stack[0]);
  }

