/* Tests the SPIntegral module. */
/* Last edited on 2005-10-02 20:52:05 by stolfi */

#define PROG_NAME "SPIntTest"

#include <SPOptions.h>
#include <SPTriang.h>
#include <SPIntegral.h>
#include <r3.h>
#include <js.h>
#include <stdio.h>

/* Computes the integral of x^n over half of the first octant. */

typedef struct Options
  { int degree;    /* Degree of polynomial to integrate. */
  } Options;

Options GetOptions(int argn, char **argc);

int main(int argn, char **argc)
  {
    Options o = GetOptions(argn, argc);
    
    auto double Func(S2Point *p);

    double Func(S2Point *p)
      { 
        int i; double r = 1.0;
        for (i = 1; i <= o.degree; i++) { r *= p->c[0]; }
        return r;
      }

    r3_t u = (r3_t){{1.0, 0.0, 0.0}};
    r3_t v = (r3_t){{0.0, 1.0, 0.0}};
    r3_t w = (r3_t){{0.0, 0.0, 1.0}};
    r3_t vw_flat, vw_sphr;

    r3_mix(0.5, &v, 0.5, &w, &vw_flat);
    r3_dir(&vw_flat, &vw_sphr);

    { /* Flat triangle integral test: */
      double sum = 0.0, corr = 0.0;
      double exact = 0.5*SQRT3/((double)(o.degree + 1)*(o.degree + 2));
      SPIntegral_OnFlatTriangle(&u, &v, &vw_flat, Func, &sum, &corr);
      fprintf(stderr, "%-15s = %+22.15e", "OnFlatTriangle", sum);
      fprintf(stderr, "  exact = %+22.15e", exact);
      fprintf(stderr, "  error = %+22.15e", sum - exact);
      fprintf(stderr, "\n");
    }
    { /* Spherical triangle integral test: */
      double sum = 0.0, corr = 0.0;
      double exact = 0.25*PI/((double)o.degree + 1);
      SPIntegral_OnSphericalTriangle(&u, &v, &vw_sphr, Func, &sum, &corr);
      fprintf(stderr, "%-15s = %+22.15e", "OnSphericalTriangle", sum);
      fprintf(stderr, "  exact = %+22.15e", exact);
      fprintf(stderr, "  error = %+22.15e", sum - exact);
      fprintf(stderr, "\n");
    }
    { /* Supersampled integral test: */
      int k; int smpOrder;
      for (k = 0, smpOrder = 1; k <= 5; k++, smpOrder *= 2)
        { S2Point_vec_t sp = (S2Point_vec_t){ 0, NULL };
          double_vec_t wp = (double_vec_t){ 0, NULL };
          int ns = 0;
          double sum = 0.0, corr = 0.0;
          double exact = 0.25*PI/((double)o.degree + 1);
          SPIntegral_SuperSampleTriangle(&u, &v, &vw_sphr, smpOrder, &sp, &wp, &ns);
          S2Point_vec_trim(&sp, ns);
          double_vec_trim(&wp, ns);
          SPIntegral_BySamples(Func, sp, wp, &sum, &corr);
          fprintf(stderr, "%-15s = %+22.15e (order = %2d)", "BySamples", sum, smpOrder);
          fprintf(stderr, "  exact = %+22.15e", exact);
          fprintf(stderr, "  error = %+22.15e", sum - exact);
          fprintf(stderr, "\n");
        }
    }
    { /* Spherical triangle area test: */
      double s = SphericalTriangleArea(&u, &v, &vw_sphr);
      double exact = 0.25*PI;
      fprintf(stderr, "%-15s = %+22.15e", "Sph.Area", s);
      fprintf(stderr, "  exact = %+22.15e", exact);
      fprintf(stderr, "  error = %+22.15e", s - exact);
      fprintf(stderr, "\n");
    }
    return 0;
  }


Options GetOptions(int argn, char **argc)
  { Options o;
    SPOptions_Parser_t *pp = SPOptions_NewParser(stderr, argn, argc);
    
    SPOptions_SetUsage(pp, 
      PROG_NAME " \\\n"
      "  -degree NUM"
    );

    SPOptions_GetKeyword(pp, "-degree");
    o.degree = SPOptions_GetNextInt(pp, 0, 10);
       
    SPOptions_Finish(pp);

    return o;
  }
