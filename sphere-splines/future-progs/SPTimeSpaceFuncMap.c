/* See SPTimeSpaceFuncMap.h. */
/* Last edited on 2005-08-28 19:12:24 by stolfi */

#include <SPBasic.h>
#include <SPTimeSpaceFuncMap.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <string.h>
#include <math.h>

#define TSFM TimeSpaceFuncMap

/* RIGHT-HAND SIDES FOR TIME-VARYING PROBLEMS: */

/* These pointwise operators are meant to be used as right-hand sides
  of the differential equation 
  
    { df(p,t) + K*(SLap(f)(p,t) - c*f(p,t)) = F(f(p,t), p, t) }
    
  where {p} is a generic point on the sphere, {t} is time, {f} is the
  unknown function, {df} is its time derivative, {SLap} is the
  spherical laplaciam operator, {K} and {c} are real coefficients, and
  {F} is a given pointwise operator.
    
  */

/*** {TSFzero(u,p,t)} *****************************************/

#define TSFzeroDescr "0"
double TSFzero(double u, S2Point *p, double t);
double TSFzero(double u, S2Point *p, double t)
  { return 0; }
static TSFM TSFzeroMap;

/*** {TSFfooo(u,p,t)} *****************************************/

#define TSFfoooDescr "3*t*x"
double TSFfooo(double u, S2Point *p, double t);
double TSFfooo(double u, S2Point *p, double t)
  { double x = p->c[0];
    return 3.0*t*x;
  }
static TSFM TSFfoooMap;

/*** {TSFexnt(u,p,t)} *****************************************/

#define TSFexntDescr "-u"
double TSFexnt(double u, S2Point *p, double t);
double TSFexnt(double u, S2Point *p, double t)
  { return (-u); }
static TSFM TSFexntMap;

/* TABLE INITIALIZATION */

static bool_t SPTimeSpaceFuncMap_initialized = FALSE;

void SPTimeSpaceFuncMap_Initialize(void);
void SPTimeSpaceFuncMap_Initialize(void)
  {
    TSFzeroMap =  (TSFM){&TSFzero, TRUE,  TSFzeroDescr};
    TSFfoooMap =  (TSFM){&TSFfooo, FALSE, TSFfoooDescr};
    TSFexntMap =  (TSFM){&TSFexnt, TRUE,  TSFexntDescr};
    SPTimeSpaceFuncMap_initialized = TRUE;
  }

SPTimeSpaceFuncMap SPTimeSpaceFuncMap_FromName(char *name)
  {
    if (! SPTimeSpaceFuncMap_initialized) { SPTimeSpaceFuncMap_Initialize(); }
    if (strcmp(name, "TSFzero") == 0) { return TSFzeroMap; }
    if (strcmp(name, "TSFfooo") == 0) { return TSFfoooMap; }
    if (strcmp(name, "TSFexnt") == 0) { return TSFexntMap; }
    affirm(FALSE, "unknown TimeSpaceFuncMap name");
    return NoTSFMap;
  }

