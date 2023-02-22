/* Test of the shcokwave method for placement optimization. */
/* Last edited on 2023-02-21 19:56:29 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include <epswr.h>
#include <bool.h>
#include <jsfile.h>
#include <jsstring.h>
#include <stmap.h>

/* INTERNAL PROTOS */

int32_t main(int32_t argc, char **argv);
void compute_forces
  ( double xini,  /* Start of first spring(fixed). */
    double *d,    /* Particle positions. */
    double *v,    /* Particle velocities. */
    double t,     /* Current time. */
    double B,     /* Internal friction coefficient. */
    double C,     /* External friction coefficient. */
    double *K,    /* Strength of each spring. */
    double *R,    /* Rest length of each spring. */
    int32_t n,    /* Number of particles. */
    double *f     /* (OUT) Force on each particle. */
  );
  
epswr_figure_t *new_figure(char *name,Interval xr, Interval yr);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    /* System parameters: */
    int32_t n = 20;      /* Number of particles. */
    double K[n];     /* Strength of each spring. */
    double R[n];     /* Rest length of each spring. */
    double m[n];     /* Particle masses. */
    double B = 16.00; /* Internal friction coeff. */
    double C = 0.00; /* External friction coeff. */

    /* System state: */
    double z[n]; /* Particle positions. */
    double v[n]; /* Particle velocities. */
    
    /* Work areas: */
    double f[n]; /* Force on each particle. */
    
    /* Position range and step: */
    double zini = 0.0;                /* Start of first spring (fixed). */
    double meanR = 1.0;               /* Mean rest length of springs. */
    double zfin = zini + 2*n*meanR;   /* Max position of last particle. */
    double zw = zfin - zini;          /* Magniude of z range. */
    
    /* Time range and step: */
    double tini = 0.0;                 /* Simulation start time. */
    double tfin = 10.0;                 /* Simulation stop time. */
    int32_t nsteps = 1000;                  /* Number of time steps. */
    double  dt = (tfin - tini)/nsteps; /* Time step. */
    
    /* Plot stuff: */
    double ts = 1.5*(zfin - zini)/(tfin - tini);   /* Plot Y coords/time. */
    Interval xr = (Interval){ zini, zfin };       /* Range of X plot coords. */
    Interval yr = (Interval){ ts*tini, ts*tfin }; /* Range of Y plot coords. */
    epswr_figure_t *eps = new_figure("out/sproing", xr, yr);
    epswr_set_pen(eps, 0.00, 0.00, 1.00,  0.15,  0, 0);

    /* System setup: */
    int32_t i;
    for(i = 0; i < n; i++) 
      { K[i] = 100.0;
        R[i] = meanR;
        m[i] = 1.00;
      }

    /* Initial state (somewhat compressed): */
    double compress = 0.5;
    for(i = 0; i < n; i++) 
      { z[i] = compress*meanR*(i+1);
        v[i] = 0;
      }
    double t = tini;
    
    /* Simulate motion: */
    int32_t it;
    for(it = 0; it < nsteps; it++)
      { /* Compute force on each particle: */
        double tnew = t + dt;
        compute_forces(zini, z, v, t, B, C, K, R, n, f);
        /* Update position and velocity, and plot: */
        for(i = 0; i < n; i++) 
          { double vnew = v[i] + dt*f[i]/m[i]; 
            double znew = z[i] + dt*(v[i] + vnew)/2.0;
            epswr_segment(eps, z[i], ts*t, znew, ts*tnew);
            epswr_set_fill_color(eps, 0.75, 0.00, 0.00);
            epswr_dot(eps, z[i], ts*t,  0.20,  TRUE, FALSE);
            if ((z[i] < zini - 2*zw) || (z[i] > zfin + 2*zw)) 
              { fprintf(stderr, "diverged\n"); it = nsteps+1; }
            z[i] = znew;
            v[i] = vnew;
           }
        t = tnew;
      }

    /* Plot final state: */
    for(i = 0; i < n; i++) 
      { epswr_set_fill_color(eps, 0.75, 0.00, 0.00); 
        epswr_dot(eps, z[i], ts*t, 0.5, TRUE, FALSE);
      }
      
    epswr_set_pen(eps, 0.00, 0.00, 0.00,  0.25,  0, 0);

    /* Add caption with parameters: */
    char buf[1024];
    epswr_text(eps, "", FALSE, 0.0, TRUE, FALSE);
    sprintf(buf, "K = %10.6f  B = %10.6f  C = %10.6f  tfin = %10.6f", K[0], B, C, tfin);
    epswr_text(eps, buf, FALSE, 0.0, TRUE, FALSE);

    epswr_frame(eps);
    epswr_end_figure(eps);
    return 0;
  }

void compute_forces
  ( double zini, 
    double *z,
    double *v,
    double t,
    double B,
    double C,
    double *K, 
    double *R, 
    int32_t n, 
    double *f
  )
  { int32_t i;
    double sfnext;  /* Force between particles {i} and {i+1}. */
    sfnext = 0.0;   /* There's no spring after the last particle. */
    for(i = n-1; i >= 0; i--)
      { /* Compute force {sf} generated by spring {i} (>0 inwards, <0 outwards): */
        double zi = z[i];
        /* Current spring length, and its time derivative: */
        double L = zi - (i > 0 ? z[i-1] : zini);
        double DL = v[i] - (i > 0 ? v[i-1] : 0.0);
        /* Relative distension (non-collapsible), and its derivative: */
        double h = 0.5*(L/R[i] - R[i]/L);
        double Dh = 0.5*DL*(1.0/R[i] + R[i]/(L*L));
        /* Force components of spring {i}: elastic and viscous: */
        double elf = K[i]*h;
        double vsf = B*Dh;
        /* Total force generated by spring {i}: */
        double sf = elf + vsf;
        /* Friction force on particle {i}: */
        double spfi = sfnext - sf;
        double frfi = -C*v[i];
        /* Total force of particle {i}: */
        f[i] = spfi + frfi;
        /* Debugging: */
        if (i == 6)
          { fprintf(stderr, "t = %8.5f\n", t); }
        if ((i == 5) || (i == 6))
          { fprintf(stderr, "  L = %10.6f  h = %10.6f  elf = %10.6f", L, h, elf);
            fprintf(stderr, "  DL = %10.6f  Dh = %10.6f  vsf = %10.6f", DL, Dh, vsf);
            fprintf(stderr, "  sf[%d] = %10.6f\n", i, sf);
          }
        if (i == 5)
          { fprintf(stderr, "  spfi = %10.6f  frfi = %10.6f\n", spfi, frfi);
            fprintf(stderr, "  f[%d] = %10.6f\n", i, f[i]);
          }
        /* Prepare for next iteration: */
        sfnext = sf;
      }
  }

epswr_figure_t *new_figure(char *name, Interval xr, Interval yr)
  {
    /* Compute ranges {xPlot, yPlot} of plot coordinates: */
    Interval xPlot = xr;  /* Range of X plot coordinate. */
    Interval yPlot = yr;  /* Range of X plot coordinate. */
    st_interval_widen(&xPlot, 1.0);
    st_interval_widen(&yPlot, 1.0);
        
    double maxSize = 150.0*epswr_pt_per_mm; /* Max figure dimensions. */
    int32_t capLines = 2;
    double capFontHeight = 10.0;
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_captioned_figure
      ( NULL, name, NULL, -1, NULL,
        xPlot.lo, xPlot.hi, yPlot.lo, yPlot.hi, maxSize, maxSize,
        capLines, capFontHeight, eps_verbose
      );
    return eps;
  }  