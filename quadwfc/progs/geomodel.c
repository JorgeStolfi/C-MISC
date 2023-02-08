/* See {geomodel.h}. */
/* Last edited on 2008-07-14 22:20:45 by stolfi */

#include <basic.h>
#include <geomodel.h>

#include <r3.h>
#include <r3x3.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>

#define GeoModel_FileFormat "2005-07-09"

void write_geomodel(FILE *wr, geomodel_t *geo)
  { auto void write_medium(medium_t *md, int i);
    auto void write_reflector(reflector_t *rf, int i);
    auto void write_bbox(interval_t bb[]);
    auto void write_interval(interval_t *r);

    int nrf = geo->nrf;
    int nmd = geo->nmd;

    int i;

    filefmt_write_header(wr, "geomodel", GeoModel_FileFormat);

    fprintf(wr, "media = %d\n", nmd);
    fprintf(wr, "reflectors = %d\n", nrf);
    fprintf(wr, "bbox = "); write_bbox(geo->bb); fprintf(wr, "\n");

    fprintf(wr, "media:\n");
    for (i = 0; i < nmd; i++)
      { 
        medium_t *md = &(geo->md[i]);
        write_medium(md, i);
      }

    fprintf(wr, "reflectors:\n");
    for (i = 0; i < nrf; i++)
      { 
        reflector_t *rf = &(geo->rf[i]);
        write_reflector(rf, i);
      }

    filefmt_write_footer(wr, "geomodel");
    fflush(wr);

    auto void write_medium(medium_t *md, int i)
      { 
        affirm(md->num == i, "invalid {num} fields");
        fprintf(wr, "%d: ", i);
        fprintf(wr, "vP = %.4f vS = %.4f", md->v_P, md->v_S);
        fputc('\n', wr);
      }
      
    auto void write_reflector(reflector_t *rf, int i)
      { 
        fprintf(wr, "%d:\n", i);
        fprintf(wr, "  mda = %d  mdb = %d\n", rf->m_above->num, rf->m_below->num);
        int nx =  rf->nv[0], ny = rf->nv[1];
        fprintf(wr, "  nx = %d ny = %d\n", nx, ny);
        fprintf(wr, "  bbox = ");
        write_bbox(rf->bb);
        fprintf(wr, "\n");
        int ix, iy;
        for (iy = 0; iy < ny; iy++)
          { 
            fprintf(wr, "  ");
            for (ix = 0; ix < nx; ix++)
              { 
                if ((ix > 0) && (ix % 5 == 0)) { fprintf(wr, "\n  "); }
                fprintf(wr, "  %.16g", rf->z[nx*iy + ix]);
              }
            fprintf(wr, "\n");
          }
      }
      
    void write_bbox(interval_t bb[])
      { 
        write_interval(&(bb[0]));
        fprintf(wr, "×"); 
        write_interval(&(bb[1]));
        fprintf(wr, "×"); 
        write_interval(&(bb[2]));
      }
    
    void write_interval(interval_t *r)
      {
        fprintf(wr, "[%.16g _ %.16g]", LO(*r), HI(*r));
      }
  }

geomodel_t read_geomodel(FILE *rd)
  { 
    auto medium_t read_medium(int i);
    auto reflector_t read_reflector(int i);
    auto void read_bbox(interval_t bb[]);
    auto interval_t read_interval(void);

    geomodel_t geo;
    filefmt_read_header(rd, "geomodel", GeoModel_FileFormat);

    int nmd = nget_int(rd, "media"); fget_eol(rd);
    affirm(nmd >= 0, "bad nmd");
    geo.nmd = nmd;
    geo.md = notnull(malloc((nmd)*sizeof(medium_t)), "no mem");

    int nrf = nget_int(rd, "reflectors"); fget_eol(rd);
    affirm(nrf >= 0, "bad nrf");
    geo.nrf = nrf;
    geo.rf = notnull(malloc((nrf+1)*sizeof(reflector_t)), "no mem");

    fget_skip_spaces(rd);
    nget_name_eq(rd, "bbox");
    read_bbox(geo.bb);
    fget_eol(rd);

    /* Read the media list: */
    fget_match(rd, "media:"); fget_eol(rd);
    int i;
    for (i = 0; i < nmd; i++)
      { 
        geo.md[i] = read_medium(i);
      }

    /* Read the reflector list: */
    fget_match(rd, "reflectors:"); fget_eol(rd);
    for (i = 0; i < nrf; i++)
      { 
        geo.rf[i] = read_reflector(i);
      }

    filefmt_read_footer(rd, "geomodel");
    return geo;

    medium_t read_medium(int i)
      { 
        medium_t md;
        md.num = fget_int(rd); 
        affirm(md.num == i, "bad medium sequence number");
        fget_skip_spaces(rd); 
        fget_match(rd, ":");
        fget_skip_spaces(rd);
        md.v_P = nget_double(rd, "vP");
        fget_skip_spaces(rd);
        md.v_S = nget_double(rd, "vS");
        fget_eol(rd);
        return md;
      }

    reflector_t read_reflector(int i)
      { 
        reflector_t rf;

        int num = fget_int(rd); 
        affirm(num == i, "bad reflector sequence number");
        fget_skip_spaces(rd); 
        fget_match(rd, ":");
        fget_eol(rd);

        fget_skip_spaces(rd);
        int mda = nget_int(rd, "mda");
        affirm((mda >= 0) && (mda < nmd), "bad medium number above");
        rf.m_above = &(geo.md[mda]);

        fget_skip_spaces(rd);
        int mdb = nget_int(rd, "mdb");
        affirm((mdb >= 0) && (mdb < nmd), "bad medium number below");
        rf.m_below = &(geo.md[mdb]);

        fget_eol(rd);

        fget_skip_spaces(rd);
        int nx = nget_int(rd, "nx");
        affirm(nx >= 0, "bad nx");
        rf.nv[0] = nx;

        fget_skip_spaces(rd);
        int ny = nget_int(rd, "ny");
        affirm(ny >= 0, "bad ny");
        rf.nv[1] = ny;

        fget_eol(rd);

        fget_skip_spaces(rd);
        nget_name_eq(rd, "bbox");
        read_bbox(rf.bb);          
        fget_eol(rd);

        rf.z = notnull(malloc(nx*ny*sizeof(double)), "no mem");
        int ix, iy;
        for (iy = 0; iy < ny; iy++)
          { for (ix = 0; ix < nx; ix++)
              { 
                fprintf(stderr, "[%d,%d]", iy, ix);
                fget_skip_formatting_chars(rd);
                rf.z[nx*iy + ix] = fget_double(rd);
              }
          }

        fget_eol(rd);
        return rf;
      }

    void read_bbox(interval_t bb[])
      { 
        int k;
        for (k = 0; k < 3; k++)
          { 
            fget_skip_spaces(rd);
            if (k > 0) { fget_match(rd, "×"); fget_skip_spaces(rd); }
            bb[k] = read_interval();
          }
      }

    interval_t read_interval(void)
      { 
        interval_t r;
        fget_match(rd, "["); 
        LO(r) = fget_double(rd);
        fget_skip_spaces(rd);
        fget_match(rd, "_"); 
        HI(r) = fget_double(rd);
        fget_skip_spaces(rd);
        fget_match(rd, "]"); 
        return r;
      }
  }

reflector_t make_sinusoidal_reflector
  (
    medium_t *m_above,
    medium_t *m_below,
    int nx, 
    int ny, 
    double xinf,
    double xsup,
    double yinf,
    double ysup,
    double zinf,
    double zsup,
    double xfreq,
    double yfreq
  )
  {
    reflector_t rf;
    /* Media: */
    rf.m_above = m_above;
    rf.m_below = m_below;

    /* Domain: */
    rf.bb[0] = (interval_t){{ xinf, xsup }}; /* X range. */
    rf.bb[1] = (interval_t){{ yinf, ysup }}; /* Y range. */
    rf.bb[2] = (interval_t){{ zinf, zsup }}; /* Z range. */
    
    /* Height field: */
    rf.nv[0] = nx;  
    rf.nv[1] = ny;
    double zctr = (zinf + zsup)/2, zamp = (zsup - zinf)/2;
    double *z = (double *)notnull(malloc(nx*ny*sizeof(double)), "no mem");
    int ix, iy;
    for (ix = 0; ix < nx; ix++)
      for (iy = 0; iy < ny; iy++)
        { 
          double xarg = ((double)ix)/((double)nx - 1) - 0.5;
          double yarg = ((double)iy)/((double)ny - 1) - 0.5;
          double arg = (2*M_PI)*(xfreq*xarg + yfreq*yarg);
          z[ix + nx*iy] = zctr + zamp*sin(arg);
        }
    rf.z = z;
    
    return rf;
  }

medium_t make_medium ( int num, double v_P, double v_S )
  {
    medium_t md;
    md.num = num;
    md.v_P = v_P;
    md.v_S = v_S;
    return md;
  }
