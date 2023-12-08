/* See {plot_utils.h} */
/* Last edited on 2023-10-01 19:38:13 by stolfi */

#include <basic.h>
#include <plot_utils.h>
#include <argparser_geo.h>

#include <r2.h>
#include <r3.h>
#include <r4.h>
#include <r3x3.h>
#include <r4x4.h>
#include <rmxn.h>
#include <hr3.h>
#include <epswr.h>
#include <affirm.h>
#include <argparser.h>
#include <bool.h>
#include <sign.h>
#include <frgb.h>

#include <values.h>
#include <stdio.h>
#include <math.h>

/* INTERNAL PROTOTYPES */

r3_t project_point(r3_t *p, hr3_pmap_t *map);
  /* Applies the perspective projection {map} to the point {p}. */
  
r3_t find_ortho_dir(r3_t *d);
  /* Finds some unit vector {u} orthogonal to {d}. */

hr3_point_t second_zenith(void);
  /* A second-best choice for the default zenith. */

bool_t is_good_zenith(hr3_point_t *upp, hr3_point_t *obs, hr3_point_t *ctr);
  /* TRUE iff {upp} is a good zenith reference for the given {obs} and {ctr}. */

frgb_t clip_color(frgb_t *c);
  /* Clips the color {c} to the unit cube, preserving 
    its brightness and hue. */

/* POSTSCRIPT FIGURE STREAMS */

epswr_figure_t *new_figure(char *name, double figSize, int32_t nCap)
  { /* Add caption only if there is a user caption. */
    /* Select a good figure size: */
    double fontSize = 14.0; /* Nominal font height (pt). */
    double mrg = 5.0;                   /* Total margin (pt). */
    /* Dimensions minus margin and captions: */ 
    double hvSize = figSize*epswr_pt_per_mm;
    /* Space for captions at bottom: */
    double capSize = nCap*fontSize; 
    double botMrg = (nCap == 0 ? mrg : mrg + capSize + mrg);
    bool_t verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, NULL, name, -1, NULL, 
        hvSize, hvSize, mrg, mrg, botMrg, mrg, verbose
      );
    epswr_set_window
      ( eps,
        mrg, mrg+hvSize, botMrg, botMrg+hvSize,
        FALSE,
        0, figSize, 0, figSize
      );
    epswr_set_text_font(eps, "Courier", fontSize);
    epswr_set_text_geometry(eps, FALSE, mrg, mrg + hvSize, mrg, mrg + capSize, 0.0);
    
    return eps;
  }
  
double default_fig_size(int32_t nRows, int32_t nCols, int32_t captionLines)
  { if (nRows <= 0) { nRows = 1; }
    if (nCols <= 0) { nCols = 1; }
    return 150.0 / (nRows < nCols ? nRows : nCols);
  }
  
/* PLOT COMPONENTS */

void draw_all_axes(epswr_figure_t *eps, hr3_pmap_t *map)
  { r3_t x = (r3_t){{1.0, 0.0, 0.0}};
    r3_t y = (r3_t){{0.0, 1.0, 0.0}};
    r3_t z = (r3_t){{0.0, 0.0, 1.0}};
    draw_axis(eps, map, &x, 1.25);
    draw_axis(eps, map, &y, 1.25);
    draw_axis(eps, map, &z, 1.25);
  }

#define AxisNBarbs 20
#define AxisHeadLength 0.150
#define AxisHeadRadius 0.015

void draw_axis
  ( epswr_figure_t *eps, 
    hr3_pmap_t *map, 
    r3_t *dir,
    double length
  )
  { 
    epswr_comment(eps, "draw_axis");
    /* Main points of arrow: */
    r3_t h, b, t; 
    r3_scale(length, dir, &h);   /* Head endpoint. */
    r3_scale(-length, dir, &t);  /* Tail endpoint. */
    r3_scale(length - AxisHeadLength, dir, &b); /* Start of arrowhead. */
    /* Project them: */
    r3_t hm = project_point(&h, map);
    r3_t tm = project_point(&t, map);
    r3_t bm = project_point(&b, map);
    /* Draw arrow shaft: */
    epswr_segment(eps, hm.c[0], hm.c[1], tm.c[0], tm.c[1]);
     /* Pick two orthogonal directions, orthogonal to {dir}: */
    r3_t u, v;
    u = find_ortho_dir(dir);
    r3_cross(dir, &u, &v);
    /* Generate arrowhead as a cone of segments: */  
    for (int32_t i = 0; i < AxisNBarbs; i++)
      { double alpha = 2*M_PI*((double)i)/((double)AxisNBarbs);
        double c = cos(alpha);
        double s = sin(alpha);
        r3_t uv;
        r3_mix(c*AxisHeadRadius, &u, s*AxisHeadRadius, &v, &uv);
        r3_t p;
        r3_add(&bm, &uv, &p);
        r3_t pm;
        pm = project_point(&p, map);
        epswr_segment(eps, tm.c[0], tm.c[1], pm.c[0], pm.c[1]);
      }
  }

r3_t project_point(r3_t *p, hr3_pmap_t *map)
  { hr3_point_t q = (hr3_point_t){{{1.0, p->c[0], p->c[1], p->c[2]}}};
    hr3_point_t r = hr3_pmap_point(&q, map);
    double w = r.c.c[0];
    return (r3_t){{r.c.c[1]/w, r.c.c[2]/w, r.c.c[3]/w}};
  }

void paint_triangle
  ( epswr_figure_t *eps,
    r3_t *P, 
    r3_t *Q, 
    r3_t *R,
    hr3_pmap_t *map,         /* Perspective projection matrix. */
    frgb_t *color,           /* Color to use. */
    r3_t *dLight,            /* Direction towards main light source. */
    double shadow            /* Amount of darkening by shadow. */
  )
  { r3_t Pm = project_point(P, map);
    r3_t Qm = project_point(Q, map);
    r3_t Rm = project_point(R, map);
    if (color->c[0] == Invisible.c[0])
      { /* invisible color */ return; }
    else
      { r3_t dMed; 
        double illum;
        frgb_t ac;
        r3_add(P, Q, &dMed);
        r3_add(R, &dMed, &dMed);
        r3_dir(&dMed, &dMed);
        illum = 1.0 + shadow * r3_cos(&dMed, dLight);
        ac = (frgb_t){{(float)(illum*color->c[0]), (float)(illum*color->c[1]), (float)(illum*color->c[2])}};
        ac = clip_color(&ac);
        epswr_set_fill_color(eps, ac.c[0], ac.c[1], ac.c[2]);
        epswr_triangle
          ( eps, 
            Pm.c[0], Pm.c[1], 
            Qm.c[0], Qm.c[1], 
            Rm.c[0], Rm.c[1], 
            TRUE, FALSE
          );
      }
  }

void draw_edge
  ( epswr_figure_t *eps,
    r3_t *P, 
    r3_t *Q, 
    hr3_pmap_t *map          /* Perspective projection matrix. */
  )
  { r3_t Pm = project_point(P, map); 
    r3_t Qm = project_point(Q, map); 
    epswr_segment(eps, Pm.c[0], Pm.c[1], Qm.c[0], Qm.c[1]);
  }

void draw_vertex
  ( epswr_figure_t *eps,
    r3_t *P, 
    hr3_pmap_t *map          /* Perspective projection matrix. */
  )
  { r3_t Pm = project_point(P, map); 
    epswr_dot(eps, Pm.c[0], Pm.c[1], 0.5, TRUE, TRUE);
  }

/* PERSPECTIVE PROJECTION */

hr3_point_t default_zenith(void)
  { return (hr3_point_t){{{0.0, 0.0, 0.0, 1.0}}}; }

hr3_point_t second_zenith(void)
  { return (hr3_point_t){{{0.0, 0.0, 1.0, 0.0}}}; }

bool_t is_good_zenith(hr3_point_t *upp, hr3_point_t *obs, hr3_point_t *ctr)
  { if ((upp == NULL) || (r4_norm(&(upp->c)) == 0.0)) { return FALSE; }
    if (upp->c.c[0] < 0.0) { return FALSE; }
    if ((obs == NULL) || (r4_norm(&(obs->c)) == 0.0))
      { return TRUE; }
    else
      { r3_t u = hr3_point_point_dir(ctr, upp);
        r3_t v = hr3_point_point_dir(ctr, obs);
        if (fabs(r3_dot(&u, &v)) >= 0.9999) { return FALSE; }
      }
    return TRUE;
  }

void fix_zenith(hr3_point_t *upp, hr3_point_t *obs, hr3_point_t *ctr, bool_t verbose)
  { hr3_point_t r;
    /* If no place to put result, do nothing: */
    if (upp == NULL) { return; }
    mumble("checking zenith reference...\n");
    if (is_good_zenith(upp, obs, ctr)) { return; }
    r = default_zenith(); 
    if (is_good_zenith(&r, obs, ctr)) { *upp = r; return; }
    r = second_zenith(); 
    if (is_good_zenith(&r, obs, ctr)) { *upp = r; return; }
    affirm(FALSE, "zenith failure");
  }

frgb_t clip_color(frgb_t *c)
  { 
    /* Get maximum coordinate {m}: */
    double m = 0.0;
    for (int32_t i = 0; i < 3; i++) { double ai = c->c[i];  if (ai > m) { m = ai; } } 
    
    /* Does it exceed the maximum? */
    if (m <= 1.0)
      { /* No -- return it: */ return *c; } 
    else
      { /* Yes -- clip it: */
        /* Find point {a} where vector {c} exits cube: */
        frgb_t a = (frgb_t){{c->c[0]/m, c->c[1]/m, c->c[2]/m}};
        /* Get brightnesses of {c} and {a}: */
        double yc = 0.3*c->c[0] + 0.6*c->c[1] + 0.1*c->c[2];
        double ya = yc/m;
        /* Add white to {a} so as to recover brightness of {c}: */
        double s = (yc < 1.0 ? yc : 1.0) - ya;
        return (frgb_t)
          {{
            a.c[0] + s*(1.0 - a.c[0]),
            a.c[1] + s*(1.0 - a.c[1]),
            a.c[2] + s*(1.0 - a.c[2])
          }};
      }
  }

/* MISCELLANEOUS */

r3_t find_ortho_dir(r3_t *d)
  { int32_t iMax, iMed, iMin;
    r3_t u;
    /* Find the max, med, and min coordinates of {d}: */
    if (fabs(d->c[0]) >= fabs(d->c[1]))
      { iMax = 0; iMed = 1; } 
    else
      { iMax = 1; iMed = 0; }
    if (fabs(d->c[2]) <= fabs(d->c[iMed]))
      { iMin = 2; }
    else if (fabs(d->c[2]) <= fabs(d->c[iMax]))
      { iMin = iMed; iMed = 2; }
    else
      { iMin = iMed; iMed = iMax; iMax = 2; }
    /* Build result: */
    u.c[iMax] = -d->c[iMed]; u.c[iMed] = d->c[iMax]; u.c[iMin] = 0.0;
    r3_dir(&u, &u);
    return u;
  }

plot_options_t parse_plot_options(argparser_t *pp)
  {
    plot_options_t o;

    /* General plotting options: */
    
    if (argparser_keyword_present(pp, "-meshSize"))
      { o.meshSize = argparser_get_next_double(pp, 0.0, 10000.0); }
    else
      { o.meshSize = 2.0; }

    o.caption = string_vec_new(10);
    { int32_t nCapLines = 0;
      while (argparser_keyword_present(pp, "-caption"))
        { string_vec_expand(&(o.caption), nCapLines);
          o.caption.e[nCapLines] = argparser_get_next(pp);
          nCapLines++;
        }
      string_vec_trim(&(o.caption), nCapLines);
    }

    if (argparser_keyword_present(pp, "-figSize"))
      { o.figSize = argparser_get_next_double(pp, 5.0, 2000.0); }
    else
      { /* Guess: user caption, plus 2 lines of figure name if not {eps}: */
        int32_t nCap = o.caption.ne;
        /* Default figure size -- fits in canvas, allowing caption: */
        o.figSize = default_fig_size(1, 1, nCap);
      }
    
    if (argparser_keyword_present(pp, "-edgeWidth"))
      { o.edgeWidth = argparser_get_next_double(pp, 0.0, DBL_MAX); }
    else
      { o.edgeWidth = 0.15; }

    /* View parameters: */
    
    if (argparser_keyword_present(pp, "-obs"))
      { o.obs.c = argparser_get_next_r4(pp, -DBL_MAX, DBL_MAX); }
    else
      { /* Invalid observer, to force computation: */
        o.obs = NoPoint;
      }

    if (argparser_keyword_present(pp, "-ctr"))
      { hr3_point_t ctrh;
        o.ctr.c = argparser_get_next_r4(pp, -DBL_MAX, DBL_MAX);
        if (o.ctr.c.c[0] <= 0.0)
          { argparser_error(pp, "center of interest must be finite"); } 
        r3_from_hr3(&ctrh);
      }
    else
      { /* Center at origin: */
        o.ctr = Origin;
      }

    if (argparser_keyword_present(pp, "-radius"))
      { o.radius = argparser_get_next_double(pp, 0.0, DBL_MAX); }
    else
      { /* Invalid radius to force computation: */
        o.radius = -1.0;
      }
    
    if (argparser_keyword_present(pp, "-up"))
      { o.upp.c = argparser_get_next_r4(pp, -DBL_MAX, DBL_MAX); }
    else
      { o.upp = Zenith; }
      
    /* Lighting parameters: */

    if (argparser_keyword_present(pp, "-light"))
      { o.light = argparser_get_next_r3_dir(pp); }
    else
      { o.light = (r3_t){{0.1, 0.2, 1.0}};
        r3_dir(&(o.light), &(o.light));
      }

    if (argparser_keyword_present(pp, "-shadow"))
      { o.shadow = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { o.shadow = 0.1; }

    return o;
  }
