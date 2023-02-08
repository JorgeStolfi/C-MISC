#ifndef irt_evray_H
#define irt_evray_H

/* Last edited on 2008-01-20 15:06:13 by stolfi */
/* Tools for evaluating a scene along a ray using IA, AA, etc. */

#include <hr3.h>
#include <r3.h>
#include <ia.h>
#include <ia_trapez.h>
#include <ia_butfly.h>
#include <aa.h>
#include <pcode.h>
#include <iaeval.h>
#include <flteval.h>

extern int irt_num_evals;  /* Total calls to the shape function  */

typedef enum 
  { arith_ia = 0,       /* Standard interval arithmentic (IA) */
    arith_ia_diff = 1,  /* Standard IA, with IA derivatives */
    arith_aa = 2,       /* Affine arithmetic (AA) */
    arith_mix = 3       /* Standard IA, or AA if result straddles 0 */
  } arith_t;
  
typedef struct shape_t  /* Description of the scene's geometry: */
  { char *proc_name;           /* Name of function that defines the scene's shape */
    pcode_proc_t proc;         /* The pseudo-code for the characteristic function */

    /* Register and stack arrays for evaluating {proc}: */
    Interval *ia_regs, *ia_stack;      /* For interval arithmetic */
    IntervalDiff *id_regs, *id_stack;  /* Ditto, including derivative. */
    AAP *aa_regs, *aa_stack;           /* For affine arithmetic. */
    Float *flt_regs, *flt_stack;       /* For floating-point. */
    FloatDiff4 *nrm_regs, *nrm_stack;  /* Ditto, including gradient. */

  } shape_t;

Float irt_eval_f_on_point
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    Float *regs,          /* Evaluation registers */
    Float *stack,         /* Evaluation stack */
    hr3_point_t *org,     /* Origin of ray */
    hr3_point_t *dst,     /* End of ray */
    Float t               /* Parameter value where to evaluate */
  );
  /* Evaluates {sh->proc} at the point of ray {org--dst} with parameter value {t},
    using ordinary floating-point. */

void irt_eval_shape_on_sub_segment
  ( shape_t *sh, 
    arith_t arith, 
    Interval *tv, 
    bool_t print_ray,
    hr3_point_t *org,
    hr3_point_t *dst,
    Interval *y, 
    ia_butfly_t *a
  );
  /* Evaluates the homogeneous four-variable function {sh} for all
    the points on the segment {org} to {dst} that correspond
    to {t} parameter values in the interval {*tv}.  

    This computation is performed with the arithmetic model
    {data->arith}.

    Returns in {y} the computed range of {f} in the interval {x}.
    Also stores in {a} a {ia_butfly_t} that contains the graph of
    {sh} along that ray segment. */
  
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
  );

void irt_seg_eval_ia_diff
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    IntervalDiff *regs,   /* Evaluation registers */
    IntervalDiff *stack,  /* Evaluation stack */
    hr3_point_t *org,     /* Start of ray */
    hr3_point_t *dst,     /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    ia_butfly_t *ip,      /* Out: bounding butterfly for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  );

void irt_seg_eval_aa
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    AAP *regs,            /* Evaluation registers */
    AAP *stack,           /* Evaluation stack */
    hr3_point_t *org,     /* Start of ray */
    hr3_point_t *dst,     /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    ia_butfly_t *ip,      /* Out: bounding butterfly for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  );

void irt_seg_eval_mix
  ( pcode_proc_t *proc,   /* The function's pseudo-code */
    Interval *iregs,      /* Evaluation registers (IA) */
    Interval *istack,     /* Evaluation stack (IA) */
    AAP *aregs,           /* Evaluation registers (AA) */
    AAP *astack,          /* Evaluation stack (AA) */
    hr3_point_t *org,     /* Start of ray */
    hr3_point_t *dst,     /* End of ray */
    Interval *tv,         /* Parameter interval (in [0__1]) */
    int print_ray,        /* Debugging flag */
    ia_butfly_t *ip,      /* Out: bounding butterfly for function along segment */
    Interval *yr          /* Out: computed range of function along segment */
  );

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
  );
  /* Divides the interval {td} into {nints} equal intervals, and
    evaluates {sh} on the corresponding segments of the {org--dst} ray,
    using {sh->arithmetic}. Returns the parameter segments in {tv[0..nints-1]},
    the function butterflys in {fvp[0..nints-1]}, and the overall
    range of the function in {fd}. */

/* NORMAL EVALUATION */

void irt_compute_surface_normal
  ( shape_t *sh,        /* The object's shape */
    arith_t arith,      /* Type of arithmetic to use */
    hr3_point_t *hit,   /* The hit point. */
    r3_t *nrm           /* Out: The normal vector at {hit}. */
  );
  /* Computes the approximate outward unit normal vector to the
    surface {sh->proc} at the point {hit}. */

/* PLOTTING */

void irt_plot_butterflies_along_segment
  ( PSStream *ps,
    int nints,           /* Number of intervals to plot */
    Interval td,         /* Overall range of parameter */
    Interval fd,         /* Overall range of function over {td} */
    Interval tv[],       /* Parameter intervals (sub-intervals of {td}) */
    ia_butfly_t fvp[] /* Function butterflies for each interval */
  );
  /* Plots a list of butterflies (such as those computed by
    {irt_eval_shape_on_equal_segments}) to the Postscript stream
    {ps}. */

void irt_plot_f_along_segment
  ( PSStream *ps,
    shape_t *sh,         /* The shape function */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *org,    /* Ray origin */
    hr3_point_t *dst,    /* Ray destination */
    int nsteps,          /* Number of intervals to plot */
    Interval td,         /* Overall range of parameter */
    Interval fd          /* Overall range of function over {td} */
  );
  /* Plots a graph of {sh->proc} along the ray {org--dst}, evaluating it
    with ordinary floating-point on {nsteps} points. */

void irt_debug_ray_graphically
  ( PSStream *ps,
    shape_t *sh,         /* The shape function */
    arith_t arith,       /* Type of arithmetic to use */
    hr3_point_t *org,    /* Ray's origin */                                  
    hr3_point_t *dst     /* Ray's destination */                             
  );
  /* Appends to {ps} a Postscript plot of the function's value along the ray
    from {org} to {dst}. */


#endif

