/* Estimate population density and natural demand at each vertex of a map. */
/* Last edited on 2023-02-21 21:36:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

#include <affirm.h>
#include <jsstring.h>
#include <jsfile.h>

#include <stmap.h>
#include <stimage.h>

/* The population density near a point {p} is estimated from the total
  length of unblocked streets in a fuzzy circle centered at {p}. The radius of
  the circle should be several city blocks.  The length of an edge is
  assumed to be its Euclidean length. */

typedef struct Options 
  { char *mapName;   /* Name of input map (minus ".rnt" extension). */
    double step;     /* Sampling step for data images (m). */
    double radius;   /* Blur radius for pop density estimate (m). */
    double totPop;   /* Total population of area covered by map. */
    char *outName;   /* Prefix for output files. */
  } Options;
  /* Command line options passed to the program. */
    
/* PROTOTYPES */

int32_t main(int32_t argc, char **argv);
Options *parse_options(int32_t argc, char **argv);
void get_next_arg_string(char **varp, int32_t *argnp, int32_t argc, char **argv, char *usage);
void get_next_arg_double(double *varp, int32_t *argnp, int32_t argc, char **argv, char *usage);
void arg_error(char *msg, char *arg, char *pname, char *usage);
Map *read_map(char *mapName);
void write_image(char *imgName, DataImage *img, double maxVal);
void write_vertex_data(char *outName, double_vec_t *d);

double_vec_t st_compute_vertex_measure(Map *m);
  /* Returns a vector {vms[0..m.nv-1]} such that {vms[vi]} is the
    street measure of vertex {vi}, defined as the total length of the
    half-edges of {m} that are incident to the vertex. */

DataImage *st_compute_street_density(Map *m, double_vec_t *vms, double step, double mrg);
  /* Computes a data image showing the approximate street density
    (meters per square meter) from a street map {m}, given the street
    measure {vms[vi]} of each vertex {vi}. The data image samples will
    be spaced {step} meters apart in each coordinate. The data image
    will cover the map's bounding box, widened by {mrg}. */

DataImage *st_estimate_pop_density(DataImage *std, double radius, double totPop);
  /* Returns a data image that gives the estimated population density
    (in people per square meter) in the neighborhood of a given point.
    The estimate is proportional to the distance-weighted average
    street density (in meters per square meter) around the point,
    where the weights are a Gaussian distribution with radius
    {radius}. The density is scaled so that the total population is
    {totPop}. */

double_vec_t st_estimate_vertex_users(Map *m, double_vec_t *vms, DataImage *pop, double totPop);
  /* Returns a vector {usr[0..m.nv-1]} such that {usr[vi]} is the
    estimated users whose nearest vertex is vertex {vi}. The users is
    assumed to be proportional to the the street measure {vms[vi]} of
    vertex {vi}, times the local population density. The proportion is
    adjusted so that the total users is {totPop}. */

int32_t main(int32_t argc, char **argv)
  { 
    Options *o = parse_options(argc, argv);
    Map *m = read_map(o->mapName);
    
    /* Compute street length measure of each vertex: */
    double_vec_t vms = st_compute_vertex_measure(m);
    write_vertex_data(txtcat(o->outName, "-vms"), &vms);
    
    /* Compute street density as function of position: */
    double mrg = 2.0*o->radius;
    DataImage *std = st_compute_street_density(m, &vms, o->step, mrg);
    Interval stdRange = st_image_sample_range(std);
    write_image(txtcat(o->outName, "-std"), std, stdRange.hi);
    
    /* Compute estimated population density as function of position: */
    fprintf(stderr, "estimating population density...\n");
    DataImage *pop = st_estimate_pop_density(std, o->radius, o->totPop);
    Interval popRange = st_image_sample_range(pop);
    fprintf(stderr, "max estimated pop density = %.2f\n", popRange.hi);
    write_image(txtcat(o->outName, "-pop"), pop, popRange.hi);
    
    /* Compute expected users at each node: */
    double_vec_t usr = st_estimate_vertex_users(m, &vms, pop, o->totPop);
    write_vertex_data(txtcat(o->outName, "-usr"), &usr);
    return 0;
  }
  
Map *read_map(char *mapName)
  { FILE *rd = open_read(txtcat(mapName, ".rnt"), TRUE);
    Map *m = st_map_read(rd);
    fclose(rd);
    return m;
  }
  
void write_image(char *imgName, DataImage *img, double maxVal)
  { FILE *wr = open_write(txtcat(imgName, ".pgm"), TRUE);
    st_image_write_as_pgm(wr, img, 0.0, maxVal);
    fclose(wr);
  }
  
void write_vertex_data(char *outName, double_vec_t *d)
  { st_write_double_vec_t(txtcat(outName, ".vdt"), "%12.6f", d); }

double_vec_t st_compute_vertex_measure(Map *m)
  { /* Assumes that street length is the Euclidean distance between endpoints: */
    /* Add street length of each vertex to each sample: */
    double_vec_t vms = double_vec_new(m->nv);
    int32_t ui;
    double totLen = 0;
    for (ui = 0; ui < m->nv; ui++)
      { VertexData *ud = m->vd[ui];
        Point *up = &(ud->p);
        /* Compute sum of length of unblocked half-edges incident to vertex {u}: */
        quad_arc_t a = m->out[ui], b = a;
        double sum = 0.0;
        do 
          { VertexData *vd = (VertexData *)quad_ddata(b);
            Point *vp = &(vd->p);
            double ruv = hypot(up->c[0] - vp->c[0], up->c[1] - vp->c[1])/2;
            sum += ruv;
            b = quad_onext(b);
          }
        while (b != a);
        /* Add the street length to the sample nearest to the vertex: */
        vms.e[ui] = sum;
        totLen += sum;
      }
    fprintf(stderr, "total street length = %.3f km\n", totLen/1000.0);
    double totArea = 50.0*totLen; /* Crude guess - improve it! */
    fprintf(stderr, "estimated urban area = %.6f km^2\n", totArea/1000000.0);
    return vms;
  }

DataImage *st_compute_street_density(Map *m, double_vec_t *vms, double step, double mrg)
  { /* Compute approximate image bounds */ 
    Interval rx, ry;
    st_map_get_bbox(m, &rx, &ry);
    st_interval_widen(&rx, mrg);
    st_interval_widen(&ry, mrg);
    DataImage *img = st_image_for_rect(rx, ry, step);
    fprintf(stderr, "data image size = %d × %d\n", img->nx, img->ny);
    /* Add street length of each vertex to each sample: */
    double S2 = step*step;  /* Area covered by each sample. */
    int32_t ui;
    for (ui = 0; ui < m->nv; ui++)
      { VertexData *ud = m->vd[ui];
        Point *up = &(ud->p);
        double msu = vms->e[ui];
        /* Add the street length to the sample nearest to the vertex: */
        st_image_increment(img, up->c[0], up->c[1], msu/S2);
      }
    return img;
  }

DataImage *st_estimate_pop_density(DataImage *std, double radius, double totPop)
  { /* Estimate relative population density from street density blurred: */
    DataImage *pop = st_image_blur(std, radius);
    /* Scale image so that total population is correct: */
    double sum = 0;
    int32_t ns = pop->nx*pop->ny, k;
    for (k = 0; k < ns; k++) { sum += pop->d[k]; }
    double step = ((double)std->mmStep)/1000.0;
    double S2 = step*step;
    double scale = totPop/(sum*S2);
    for (k = 0; k < ns; k++) { pop->d[k] *= scale; }
    return pop;
  }

double_vec_t st_estimate_vertex_users(Map *m, double_vec_t *vms, DataImage *pop, double totPop)
  { double_vec_t usr = double_vec_new(m->nv);
    /* Compute estimated users per node, unscaled. */
    int32_t ui;
    double sum = 0.0;
    for (ui = 0; ui < m->nv; ui++)
      { VertexData *ud = m->vd[ui];
        Point *up = &(ud->p);
        double popu = st_image_interpolate(pop, up->c[0], up->c[1]);
        usr.e[ui] = popu*vms->e[ui];
        sum += usr.e[ui];
      }
    /* Scale users vector to give the specified total users. */
    double scale = totPop/sum;
    for (ui = 0; ui < m->nv; ui++) { usr.e[ui] *= scale; }
    return usr;
  }
  
#define ARG_ERROR(Msg,Arg) arg_error((Msg),(Arg),argv[0],usage)
  
#define GET_NEXT_STRING(Var) get_next_arg_string(&(Var), &argn, argc, argv, usage)
#define GET_NEXT_DOUBLE(Var) get_next_arg_double(&(Var), &argn, argc, argv, usage)
  
Options *parse_options(int32_t argc, char **argv)
   /*
     Parses the command line options, returns a record with their options. */
  {
    char* usage = "\n  [ -help ] -mapName NAME \\\n  -step NUM -radius NUM -totPop NUM \\\n  -outName NAME";
    Options *o = (Options *)malloc(sizeof(Options));
    int32_t argn;

    /* Defaults: */
    o->mapName = NULL; /* For now. */
    o->step = -1.0;    /* For now. */
    o->radius = -1.0;  /* For now. */
    o->totPop = -1.0;  /* For now. */
    o->outName = NULL; /* For now. */
    
    argn = 1;

    /* Scan command line options. */
    while ((argn < argc) && (argv[argn][0] == '-') && (argv[argn][1] != '\0'))
      {
        char *key = argv[argn];
        if ((key[0] == '-') && (key[1] == '-') && (key[2] != '\0')) { key++; }
        if (strcmp(key, "-help") == 0)
          { fprintf(stderr, "usage: %s %s\n", argv[0], usage); exit(0); }
        else if (strcmp(key, "-mapName") == 0)
          { GET_NEXT_STRING(o->mapName);
          }
        else if (strcmp(key, "-step") == 0)
          { GET_NEXT_DOUBLE(o->step);
          }
        else if (strcmp(key, "-radius") == 0)
          { GET_NEXT_DOUBLE(o->radius);
          }
        else if (strcmp(key, "-totPop") == 0)
          { GET_NEXT_DOUBLE(o->totPop);
          }
        else if (strcmp(key, "-outName") == 0)
          { GET_NEXT_STRING(o->outName);
          }
        else 
          { ARG_ERROR("unknown option", argv[argn]); }
        ++argn;
      }
    if (argn != argc) { ARG_ERROR("extraneous arguments", argv[argn]); }
    if (o->mapName == NULL) { ARG_ERROR("must specify \"-mapName\"", ""); }
    if (o->step < 0.0) { ARG_ERROR("must specify \"-step\"", ""); }
    if (o->radius < 0.0) { ARG_ERROR("must specify \"-radius\"", ""); }
    if (o->totPop < 0.0) { ARG_ERROR("must specify \"-totPop\"", ""); }
    if (o->outName == NULL) { ARG_ERROR("must specify \"-outName\"", ""); }
    return o;
  }

void get_next_arg_string(char **varp, int32_t *argnp, int32_t argc, char **argv, char *usage)
   /*
     Stores the next command line argument (as a string) into "*varp" */
  {
    int32_t argn = *argnp;
    if (argn+1 >= argc)
      { ARG_ERROR("missing arg value", argv[argn]); }
    (*varp) = argv[argn+1];
    (*argnp) = argn+1;
  }
  
void get_next_arg_double(double *varp, int32_t *argnp, int32_t argc, char **argv, char *usage)
   /*
     Stores the next command line argument (as a double) into "*varp" */
  {
    int32_t argn = *argnp;
    char *end;
    if (argn+1 >= argc)
      { ARG_ERROR("missing arg value", argv[argn]); }
    (*varp) = strtod(argv[argn+1], &end);
    if ((*end) != '\0') 
      { ARG_ERROR("invalid numeric argument", argv[argn+1]); }
    (*argnp) = argn+1;
  }

void arg_error(char *msg, char *arg, char *pname, char *usage)
   /*
     Prints "msg", "arg", the "usage" message, and exits.
     Handy while parsing the command line arguments. */
  {
    fprintf(stderr, "%s %s\n", msg, arg);
    fprintf(stderr, "usage: %s %s\n", pname, usage);
    exit(1);
  }
