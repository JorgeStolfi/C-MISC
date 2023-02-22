/* Find optimum placement of public telephones in a street map. */
/* Last edited on 2023-02-22 09:28:48 by stolfi */

#include <stmap.h>
#include <stimage.h>
#include <tup.h>
#include <affirm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <quad.h>
#include <intsort.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>

/* Optimize the placement of public telephones in a street map.
  Assumes that the edge costs are symmetric.
*/

typedef struct Options 
  { char *mapName;        /* Name of input map (minus ".rnt" extension). */
    char *demName;        /* Name of basic demand file (minus ".vdt" extension). */
    double maxDist;       /* Path-cost to nearest phone must not exceed this. */
    double partDist;      /* Path-cost at which half of the demand is lost. */
    double maxCapture;    /* Demand capture fraction for zero path-cost. */
    char *outName;        /* Prefix for output files. */
    /* Parameters for strippy heuristic: */
    int32_t nStripWidths;     /* Number of strip widthds to try. */
    double minStripWidth; /* Min strip width to try. */
    double maxStripWidth; /* Min strip width to try. */
  } Options;
  /* Command line options passed to the program. */

/* PROTOTYPES */

int32_t main(int32_t argc, char **argv);
Options *parse_options(int32_t argc, char **argv);
void get_next_arg_string(char **varp, int32_t *argnp, int32_t argc, char **argv, char *usage);
void get_next_arg_double(double *varp, int32_t *argnp, int32_t argc, char **argv, char *usage);
void get_next_arg_int32_t(int32_t *varp, int32_t *argnp, int32_t argc, char **argv, char *usage);
void arg_error(char *msg, char *arg, char *pname, char *usage);
Map *read_map(char *mapName);
double_vec_t read_vertex_data(char *demName);
void write_phones(char *phoName, Map *m, phone_vec_t *ph);

phone_vec_t st_honey_place_phones(Map *m, double maxDist);
  /* Attempts to place phones on the street map {m} so that 
    every vertex is {maxDist}-reachable from some phone.
    Uses the hexagonal (honeycomb) heuristics.  May leave 
    parts of the map uncovered. */

void st_greedy_complete_phones(Map *m, double maxDist, phone_vec_t *ph);
  /* Appends to the vector {*ph} zero or more new phones, so that all
    vertices of {m} are {maxDist}-reachable from some phone.
    Uses a greedy heuristic, processing the vertices in a fixed order. */
    
void st_strippy_complete_phones(Map *m, double maxDist, double stripWd, phone_vec_t *ph);
  /* Appends to the vector {*ph} zero or more new phones, so that all
    vertices of {m} are {maxDist}-reachable from some phone.
    Uses a strip-coverage heuristic, with vertical strips of width {stripWd}. */

int32_t st_greedy_cover_vertex(Map *m, int32_t vi, double maxDist, int32_t *vcover, int32_t *ecover);
  /* Returns a vertex {wi} that is {maxDist}-reachable from {vi} and
    is such that placing an additional phone at {wi} will maximizes
    the number of vertices that are {maxDist}-reachable from some
    phone. Also increments {vcover[ui]} and {ecover[ei]} for all
    vertices {vi} and edges {ei} in the {maxDist}-ball of {wi}. */

int32_t st_cover_vertices_in_strip
  ( Map *m,           /* Map. */
    int32_t *vseq,      /* Unprocessed vertices in strip, in Y order. */
    int32_t nvseq,      /* Number of vertices in {vrem}. */ 
    double maxDist, /* Maxmum path-cost. */ 
    int32_t *vcover,    /* Vertex coverage by {maxDist}-balls of fixed phones. */ 
    int32_t *ecover     /* Edge coverage by {maxDist}-balls of fixed phones. */ 
  );
  /* Returns a vertex {wi} such that placing one additional phone at
    {wi} will maximize the number {nc} of consecutive vertices
    {vseq[0..nc]} that are {maxDist}-reachable from some phone. Also
    increments {vcover[ui]} and {ecover[ei]} for all vertices {vi} and
    edges {ei} in the {maxDist}-ball of {wi}. */
  
int32_t main(int32_t argc, char **argv)
  { Options *o = parse_options(argc, argv);
    /* Read data: */
    Map *m = read_map(o->mapName);
    double_vec_t dem = read_vertex_data(o->demName);
    affirm(dem.ne == m->nv, "inconsistent vertex count");
    
    /* The phone placement: */
    phone_vec_t ph;   /* phone_t layout. */
    
    /* Lost demand per vertex: */
    double_vec_t lost = double_vec_new(m->nv);
    
    /* Compute strip placement for various widths, select best: */
    int32_t nwd = o->nStripWidths;
    double magmin = o->minStripWidth;
    double magmax = o->maxStripWidth;
    double bestwd = 0.0;
    int32_t bestnph = 1000000;
    int32_t t;
    for (t = 0; t <= nwd + 1; t++)
      { double tf = ((double)t)/((double)nwd);
        double stripWd;
        if (t <= nwd) 
          { double mag = magmin*exp(tf*log(magmax/magmin));
            stripWd = mag * o->maxDist;
          }
        else
          { stripWd = bestwd; }
        ph = phone_vec_new(0);
        st_strippy_complete_phones(m, o->maxDist, stripWd, &ph);
        fprintf
          ( stderr, "strippy heuristic (mag = %5.3f): %d phones\n",
            stripWd/o->maxDist, ph.ne
          );
        if (t <= nwd)
          { if (ph.ne < bestnph) 
              { bestnph = ph.ne; bestwd = stripWd; }
          }
        else 
          { st_recompute_phone_usage(m, &dem, o->partDist, o->maxCapture, &ph, &lost);
            write_phones(txtcat(o->outName, "-strip"), m, &ph);
            free(ph.e);
          }
      }
    
    /* Compute greedy placement: */
    ph = phone_vec_new(0);
    st_greedy_complete_phones(m, o->maxDist, &ph);
    fprintf(stderr, "greedy heuristic: %d phones\n", ph.ne);
    st_recompute_phone_usage(m, &dem, o->partDist, o->maxCapture, &ph, &lost);
    write_phones(txtcat(o->outName, "-greed"), m, &ph);
    free(ph.e);
    
    /* Compute honeycomb placement: */
    ph = st_honey_place_phones(m, o->maxDist);
    fprintf(stderr, "hexag heuristic: %d phones\n", ph.ne);
    /* Fill holes with strippy heuristic: */
    { double stripWd = o->maxDist*sqrt(2.0);
      st_strippy_complete_phones(m, o->maxDist, stripWd, &ph);
      fprintf(stderr, "hexag+strippy heuristic: %d phones\n", ph.ne);
    }
    st_recompute_phone_usage(m, &dem, o->partDist, o->maxCapture, &ph, &lost);
    write_phones(txtcat(o->outName, "-hexag"), m, &ph);
    free(ph.e);
    
    return 0;
  }
  
phone_vec_t st_honey_place_phones(Map *m, double maxDist)
  { 
    /* Ideally, the positions of the phones is such that an array of disks
    with the specified {radius}, centered at those points, will cover
    the bounding box of {m}, with zero slop.
    
    Each phone is actually assigned to the vertex of {m} that is
    closest to its ideal position above (in the Euclidean metric),
    provided that the distance between the two points does not exceed
    {pMax}. If the distance is greater than {pMax}, that phone is
    omitted. */
    
    /* Preserve old phones: */
    phone_vec_t ph = phone_vec_new(10);
    int32_t nph = 0;        /* Phones in layout are {ph[0..nph-1]}. */
    /* Compute map bounds */ 
    Interval rx, ry;
    st_map_get_bbox(m, &rx, &ry);
    /* Maximum allowed deviation of phone vertex from ideal phone position: */ 
    double maxJitter = 51.0;       
    /* Most-cases maximum deviation of phone vertex from ideal phone position: */ 
    double expJitter = 10.0;       
    /* Radius of ideal cover disks (centered on ideal positions): */ 
    double hRadius = maxDist/sqrt(2.0) - expJitter;
    /* Displacement between phones on same honeycomb row: */ 
    double xStep = sqrt(3.0) * hRadius; 
    /* Displacement between honeycomb rows: */ 
    double yStep = 1.5 * hRadius; 
    /* Number of honeycomb rows: */
    int32_t nrows = (int32_t)ceil((ry.hi - ry.lo - hRadius)/yStep) + 1;
    /* Overcoverage of Y range: */
    double ySlop = yStep*(nrows-1) + hRadius - (ry.hi - ry.lo);
    /* Y position of first row: */
    double y = ry.lo + hRadius/2 - ySlop/2;
    /* Cover rectangle {rx × ry} with honeycomb of phones: */
    int32_t row = 0;
    while(y - hRadius/2 < ry.hi)
      { double x = rx.lo + (row % 2 == 0 ? 0.0 : xStep/2);
        while (x - xStep/2 < rx.hi)
          { Point ctr = (Point){{x, y}}; /* Ideal phone location. */
            int32_t vi = st_map_nearest_vertex(m, ctr); /* Vertex assigned to phone. */
            VertexData *vd = m->vd[vi];
            Point loc = vd->p;           /* Assigned position of phone. */
            double disp = hypot(ctr.c[0] - loc.c[0], ctr.c[1] - loc.c[1]);
            if (disp <= maxJitter) 
              { phone_vec_expand(&ph, nph);
                phone_t *phi = &(ph.e[nph]);
                phi->vertex = vi;
                phi->usage = -1.0; /* Unknown. */
                nph++;
              }
            x += xStep;
          }
        row++; 
        y += yStep;
      }
    phone_vec_trim(&ph, nph);
    return ph;
  }
  
void st_greedy_complete_phones(Map *m, double maxDist, phone_vec_t *ph)
  { 
    int32_t nph = ph->ne;
    /* Compute coverage of existing phones: */
    int32_t vcover[m->nv];    /* Vertex coverage by TUP path-cost balls. */
    int32_t ecover[m->ne];    /* Edge coverage by TUP path-cost balls. */
    int32_t uph[nph];
    float dMax[nph];
    int32_t i;
    for (i = 0; i < nph; i++)
      { /* Find vertex nearest to specified center: */
        uph[i] = ph->e[i].vertex;
        dMax[i] = (float)maxDist; 
      }
    st_compute_coverage(m, uph, dMax, nph, vcover, ecover);

    /* Make list of vertices sorted by increasing X coord: */

    auto int32_t xyOrder(int32_t ui, int32_t vi);
      /* -1 iff vertex {ui} precedes {vi} in X-then-Y lex order. */
    
    int32_t xyOrder(int32_t ui, int32_t vi)
      { VertexData *ud = m->vd[ui];
        VertexData *vd = m->vd[vi];
        double ux = ud->p.c[0], uy = ud->p.c[1];
        double vx = vd->p.c[0], vy = vd->p.c[1];
        /* Compare X coords: */
        if (ux < vx)
          { return -1; }
        else if (ux > vx)
          { return +1; }
        /* Compare Y coords: */
        if (uy < vy)
          { return -1; }
        else if (uy > vy)
          { return +1; }
        /* Compare vertex indices: */
        if (ui < vi)
          { return -1; }
        else if (ui > vi)
          { return +1; }
        /* Same vertex! */
        return 0;
      }

    int32_t vseq[m->nv];
    int32_t k;
    for (k = 0; k < m->nv; k++) { vseq[k] = k; }
    isrt_heapsort(vseq, m->nv, xyOrder, +1);
    /* Cover uncovered vertices: */
    for (k = 0; k < m->nv; k++)
      { int32_t vi = vseq[k];
        if (vcover[vi] == 0)
          { int32_t wi = st_greedy_cover_vertex(m, vi, maxDist, vcover, ecover); 
            phone_vec_expand(ph, nph);
            ph->e[nph] = (phone_t){ wi, -1.0 };
            nph++;
          }
      }
    phone_vec_trim(ph, nph);
  }

void st_strippy_complete_phones(Map *m, double maxDist, double stripWd, phone_vec_t *ph)
  { 
    /* The map is divided into `strips' parallel to the Y axis, whose
      width is {stripWd}.  Within each strip, vertices are covered from 
      low Y to high Y. */
    
    /* Count of installed phones: */
    int32_t nph = ph->ne;

    /* Compute map bounds */ 
    Interval rx, ry;
    st_map_get_bbox(m, &rx, &ry);
    /* Number of strips needed to cover the map: */
    int32_t nstrips = (int32_t)floor((rx.hi - rx.lo + 0.00001)/stripWd) + 1;
    /* Overcoverage of X range: */
    double xSlop = stripWd*nstrips - (ry.hi - ry.lo);
    /* Low X position of first strip: */
    double xStart = rx.lo - xSlop/2;
    
    /* Compute strip index {strip[ui]} of each vertex: */
    int32_t strip[m->nv];
    int32_t ui;
    for (ui = 0; ui < m->nv; ui++)
      { VertexData *ud = m->vd[ui];
        double ux = ud->p.c[0];
        int32_t us = (int32_t)floor((ux - xStart)/stripWd);
        if (us < 0) { us = 0; }
        if (us >= nstrips) { us = nstrips - 1; }
        strip[ui] = us;
      }
    
    /* Compute coverage of existing phones: */
    int32_t vcover[m->nv];    /* Vertex coverage by TUP path-cost balls. */
    int32_t ecover[m->ne];    /* Edge coverage by TUP path-cost balls. */
    { int32_t uph[nph];
      float dMax[nph];
      int32_t i;
      for (i = 0; i < nph; i++)
        { /* Find vertex nearest to specified center: */
          uph[i] = ph->e[i].vertex;
          dMax[i] = (float)maxDist; 
        }
      st_compute_coverage(m, uph, dMax, nph, vcover, ecover);
    }
    
    /* Make list of vertices sorted by strip then by increasing Y: */
    
    auto int32_t stripOrder(int32_t ui, int32_t vi);
      /* -1 iff vertex {ui} lies on an earlier strip than {vi}. */
    
    int32_t stripOrder(int32_t ui, int32_t vi)
      { VertexData *ud = m->vd[ui];
        int32_t us = strip[ui];
        double uy = ud->p.c[1];
        
        VertexData *vd = m->vd[vi];
        int32_t vs = strip[vi];
        double vy = vd->p.c[1];

        /* Compare strip indices: */
        if (us < vs)
          { return -1; }
        else if (us > vs)
          { return +1; }
        /* Compare Y coords: */
        if (uy < vy)
          { return -1; }
        else if (uy > vy)
          { return +1; }
        /* Compare vertex indices: */
        if (ui < vi)
          { return -1; }
        else if (ui > vi)
          { return +1; }
        /* Same vertex! */
        return 0;
      }
    
    int32_t vseq[m->nv];
    int32_t k;
    for (k = 0; k < m->nv; k++) { vseq[k] = k; }
    isrt_heapsort(vseq, m->nv, stripOrder, +1);

    /* Conver uncovered vertices in strip order: */
    int32_t ini = 0; /* Next unprocessed vertex is {vseq[ini]}. */
    while(ini < m->nv)
      { /* First vertex of current strip is {vseq[ini]}. */
        int32_t s = strip[vseq[ini]];
        /* Count of vertices in current `bite' of strip. */
        int32_t nbite = 0; 
        /* Find first vertex {vseq[lim]} of next strip: */
        int32_t lim = ini + 1;
        while ((lim < m->nv) && (strip[vseq[lim]] == s)) { lim++; }
        affirm((lim >= m->nv) || (strip[vseq[lim]] > s), "strip sort bug");
        /* Try to cover strip in Y order: */
        while (ini < lim)
          { int32_t vi = vseq[ini];
            if (vcover[vi] == 0)
              { /* End previous bite: */
                if (nbite > 0) { fprintf(stderr, "[%d]", nbite); }
                /* Cover as many consec. vtces. as possible starting with {vi}: */
                int32_t *vrem = &(vseq[ini]); /* Unprocessed vertices in strip. */
                int32_t nvrem = lim - ini;    /* Number of such vertices. */
                int32_t wi = st_cover_vertices_in_strip
                  ( m, vrem, nvrem, maxDist, vcover, ecover ); 
                affirm(vcover[vi] > 0, "{vi} not covered");
                phone_vec_expand(ph, nph);
                ph->e[nph] = (phone_t){ wi, -1.0 };
                nph++;
                /* Start new bite: */ 
                nbite = 1;                
              }
            else
              { nbite++; }
            ini++;
          }
        /* End last bite: */
        if (nbite > 0) { fprintf(stderr, "[%d]\n", nbite); }
        /* Prepare for next strip: */
        ini = lim;
      }
    phone_vec_trim(ph, nph);
  }

int32_t st_greedy_cover_vertex(Map *m, int32_t vi, double maxDist, int32_t *vcover, int32_t *ecover)
  { 
    /* Work files for {st_map_compute_costs}: */
    float d[m->nv];
    float c[2 * m->ne];
    quad_arc_t e[m->nv];
    float dMax = (float)maxDist;
    /* Candidates for TUP placement: */
    int32_t cand[m->nv];
    int32_t ncand;
    /* The following
      implementation assumes that the graph {m} is symmetric.
      Otherwise we should use reverse edges when obtaining the
      candidate list {cand}. */
    /* The candidates are all vertices {dMax}-reachable from {vi}: */ 
    st_map_init_costs(m, d, e, c);
    st_map_compute_costs(m, vi, dMax, cand, &ncand, d, e, c);
    /* Reset tables for coverage computations: */
    st_map_reset_costs(m, cand, ncand, d, e, c);
    /* Among those candidates, select the one with maximum `new' coverage: */
    int32_t k;
    int32_t bi = -1, bnc = -1;
    for (k = 0; k < ncand; k++)   
      { int32_t wi = cand[k];
        /* Find coverage of {wi}: */
        int32_t r[m->nv];
        int32_t nr;
        st_map_compute_costs(m, wi, dMax, r, &nr, d, e, c);
        /* Count number of vertices that are excluvely covered by {wi}: */
        int32_t nc = 0, j;
        for (j = 0; j < nr; j++)
          { int32_t ui = r[j];
            if (vcover[ui] == 0) { nc++; }
          }
        /* Remember best candidate: */
        if (nc > bnc) { bnc = nc; bi = wi; }
        /* Reset tables for next candidate: */
        st_map_reset_costs(m, r, nr, d, e, c);
      }
    /* Update {vcover,ecover} for best canadidate: */
    { /* Find {dMax}-ball of {bi}: */
      int32_t r[m->nv];
      int32_t nr;
      st_map_compute_costs(m, bi, dMax, r, &nr, d, e, c);
      /* Increment coverage of elems in ball: */
      st_increment_coverage(m, dMax, r, nr, d, c, vcover, ecover);
      /* Reset tables, just to keep the habit: */
      st_map_reset_costs(m, r, nr, d, e, c);
    }
    return bi;
  }

int32_t st_cover_vertices_in_strip
  ( Map *m,           /* Map. */
    int32_t *vseq,      /* Unprocessed vertices in strip, in Y order. */
    int32_t nvseq,      /* Number of vertices in {vrem}. */ 
    double maxDist, /* Maxmum path-cost. */ 
    int32_t *vcover,    /* Vertex coverage by {maxDist}-balls of fixed phones. */ 
    int32_t *ecover     /* Edge coverage by {maxDist}-balls of fixed phones. */ 
  )
  { affirm(nvseq > 0, "no vertex to cover??"); 
    /* Vertex that MUST be covered: */
    int32_t vi = vseq[0];
    affirm(vcover[vi] == 0, "{vi} already covered"); 
    /* Work files for {st_map_compute_costs}: */
    float d[m->nv];
    float c[2 * m->ne];
    quad_arc_t e[m->nv];
    float dMax = (float)maxDist;
    /* Candidates for TUP placement: */
    int32_t cand[m->nv];
    int32_t ncand;
    /* The following
      implementation assumes that the graph {m} is symmetric.
      Otherwise we should use reverse edges when obtaining the
      candidate list {cand}. */
    /* The candidates are all vertices {dMax}-reachable from {vi}: */ 
    st_map_init_costs(m, d, e, c);
    st_map_compute_costs(m, vi, dMax, cand, &ncand, d, e, c);
    /* Reset tables for coverage computations: */
    st_map_reset_costs(m, cand, ncand, d, e, c);
    /* Among those candidates, select the best one: */
    int32_t candk;
    int32_t bi = -1, bscore = -1, buntie = -1;
    for (candk = 0; candk < ncand; candk++)   
      { int32_t wi = cand[candk]; /* Next candidate for phone location. */
        /* Find {maxDist}-ball of {wi}: */
        int32_t r[m->nv];
        int32_t nr;
        st_map_compute_costs(m, wi, dMax, r, &nr, d, e, c);
        /* Increment coverage of vertices in {maxDist}-ball of {wi}: */
        int32_t untie = 0; /* Tie-breaking score: number of newly covered vtces. */
        { int32_t j; 
          for (j = 0; j < nr; j++) 
            { vcover[r[j]]++; if (vcover[r[j]] == 1) { untie++; } }
        }
        /* Compute candidate's {score} and tie-breaking criterion: */
        int32_t score = 0; /* Score of {wi}: index of first uncovered vtx. */
        while((score < nvseq) && (vcover[vseq[score]] > 0)) { score++; }
        affirm(vcover[vi] > 0, "{vi} not covered");
        affirm(score > 0, "scoring bug");
        /* Remember best candidate: */
        if ((score > bscore) || ((score == bscore) && (untie > buntie)))
          { bi = wi; bscore = score; buntie = untie; }
        /* Restore vertex coverage counts: */
        { int32_t j; for (j = 0; j < nr; j++) { vcover[r[j]]--; } }
        /* Reset distance tables for next candidate: */
        st_map_reset_costs(m, r, nr, d, e, c);
      }
    /* Update {vcover,ecover} for best canadidate: */
    { /* Find {dMax}-ball of {bi}: */
      int32_t r[m->nv];
      int32_t nr;
      st_map_compute_costs(m, bi, dMax, r, &nr, d, e, c);
      /* Increment coverage of elems in ball: */
      st_increment_coverage(m, dMax, r, nr, d, c, vcover, ecover);
      /* Reset tables, just to keep the habit: */
      st_map_reset_costs(m, r, nr, d, e, c);
    }
    return bi;
  }

Map *read_map(char *mapName)
  { FILE *rd = open_read(txtcat(mapName, ".rnt"), TRUE);
    Map *m = st_map_read(rd);
    fclose(rd);
    return m;
  }
  
void write_phones(char *phoName, Map *m, phone_vec_t *ph)
  { st_write_phones(txtcat(phoName, ".phn"), m, ph); } 

double_vec_t read_vertex_data(char *demName)
  { return st_read_double_vec_t(txtcat(demName, ".vdt")); }
  
#define ARG_ERROR(Msg,Arg) arg_error((Msg),(Arg),argv[0],usage)
  
#define GET_NEXT_STRING(Var) get_next_arg_string(&(Var), &argn, argc, argv, usage)
#define GET_NEXT_DOUBLE(Var) get_next_arg_double(&(Var), &argn, argc, argv, usage)
#define GET_NEXT_INT32_T(Var) get_next_arg_int32_t(&(Var), &argn, argc, argv, usage)
  
Options *parse_options(int32_t argc, char **argv)
   /*
     Parses the command line options, returns a record with their options. */
  {
    char* usage = "\n  [ -help ] -mapName NAME \\\n  -maxDist NUM -partDist NUM -maxCapture NUM \\\n  -outName NAME";
    Options *o = (Options *)malloc(sizeof(Options));
    int32_t argn;

    /* Defaults: */
    o->mapName = NULL;      /* For now. */
    o->demName = NULL;      /* For now. */
    o->maxDist = -1.0;      /* For now. */
    o->partDist = -1.0;     /* For now. */
    o->maxCapture = -1.0;   /* For now. */
    o->nStripWidths = -1;   /* For now. */
    o->outName = NULL;      /* For now. */
    
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
        else if (strcmp(key, "-demName") == 0)
          { GET_NEXT_STRING(o->demName);
          }
        else if (strcmp(key, "-maxDist") == 0)
          { GET_NEXT_DOUBLE(o->maxDist);
            if (o->maxDist <= 0.0)
              { ARG_ERROR("invalid \"-maxDist\"", ""); }
          }
        else if (strcmp(key, "-partDist") == 0)
          { GET_NEXT_DOUBLE(o->partDist);
            if (o->partDist <= 0.0)
              { ARG_ERROR("invalid \"-partDist\"", ""); }
          }
        else if (strcmp(key, "-maxCapture") == 0)
          { GET_NEXT_DOUBLE(o->maxCapture);
            if ((o->maxCapture <= 0.0) || (o->maxCapture > 0.999))
              { ARG_ERROR("invalid \"-maxCapture\"", ""); }
          }
        else if (strcmp(key, "-outName") == 0)
          { GET_NEXT_STRING(o->outName);
          }
        else if (strcmp(key, "-stripWidth") == 0)
          { double stripWidth;
            GET_NEXT_DOUBLE(stripWidth);
            if ((stripWidth <= 0.0) || (stripWidth > 10.00))
              { ARG_ERROR("invalid \"-stripWidth\"", ""); }
            o->minStripWidth = stripWidth;
            o->maxStripWidth = stripWidth;
            o->nStripWidths = 1;
          }
        else if (strcmp(key, "-stripWidths") == 0)
          { GET_NEXT_INT32_T(o->nStripWidths);
            if ((o->nStripWidths <= 0) || (o->nStripWidths > 100.0))
              { ARG_ERROR("invalid \"-stripWidths\" (count)", ""); }
            GET_NEXT_DOUBLE(o->minStripWidth);
            if ((o->minStripWidth <= 0.0) || (o->minStripWidth > 10.00))
              { ARG_ERROR("invalid \"-stripWidths\" (min)", ""); }
            GET_NEXT_DOUBLE(o->maxStripWidth);
            if ((o->maxStripWidth <= 0.0) || (o->maxStripWidth > 10.00))
              { ARG_ERROR("invalid \"-stripWidths\" (max)", ""); }
          }
        else 
          { ARG_ERROR("unknown option", argv[argn]); }
        ++argn;
      }
    if (argn != argc) { ARG_ERROR("extraneous arguments", argv[argn]); }
    if (o->mapName == NULL) { ARG_ERROR("must specify \"-mapName\"", ""); }
    if (o->demName == NULL) { ARG_ERROR("must specify \"-demName\"", ""); }
    if (o->maxDist < 0.0) { ARG_ERROR("must specify \"-maxDist\"", ""); }
    if (o->partDist < 0.0) { ARG_ERROR("must specify \"-partDist\"", ""); }
    if (o->maxCapture < 0.0) { o->maxCapture = 0.95; }
    if (o->outName == NULL) { ARG_ERROR("must specify \"-outName\"", ""); }
    if (o->nStripWidths < 0) 
      { o->nStripWidths = 1;
        o->minStripWidth = 0.50;
        o->maxStripWidth = 0.50;
      }
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

void get_next_arg_int32_t(int32_t *varp, int32_t *argnp, int32_t argc, char **argv, char *usage)
   /*
     Stores the next command line argument (as an int32_t) into "*varp" */
  {
    int32_t argn = *argnp;
    double v;
    char *end;
    if (argn+1 >= argc)
      { ARG_ERROR("missing arg value", argv[argn]); }
    v = strtod(argv[argn+1], &end);
    if ((*end) != '\0') 
      { ARG_ERROR("invalid numeric argument", argv[argn+1]); }
    (*varp) = (int32_t)(v + 0.5);
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

