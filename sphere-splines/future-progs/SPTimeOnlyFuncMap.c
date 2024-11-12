/* See SPTimeOnlyFuncMap.h. */
/* Last edited on 2005-08-21 15:10:10 by stolfi */

#include <SPBasic.h>
#include <SPTimeOnlyFuncMap.h>
#include <affirm.h>
#include <nat.h>
#include <bool.h>
#include <string.h>
#include <math.h>

#define TOFM TimeOnlyFuncMap

/* RIGHT-HAND SIDES FOR TIME-VARYING PROBLEMS: */

/* These pointwise operators are meant to be used as right-hand sides
  of the differential equation 
  
    { M·TDD(f)(t) + R·TD(f)(t) + K·f(t) = TMap(f(t), t) }
    
  where {f} is the function to be determined, which depends on time
  {t}; {TD} denotes single derivative with respect to time;
  {TDD} denotes double derivative with respect to time;
  {M}, {R}, {K} are constants.    
*/

#define TFzeroDescr "0"
double TFzero(double u, double t);
double TFzero(double u, double t)
  { return 0.0; }
static TOFM TFzeroMap;

#define TFstepDescr "(t>0 ? 1 : 0)"
double TFstep(double u, double t);
double TFstep(double u, double t)
  { return (t > 0 ? 1.0 : 0.0); }
static TOFM TFstepMap;

#define TFselfDescr "u"
double TFself(double u, double t);
double TFself(double u, double t)
  { return u; }
static TOFM TFselfMap;

#define TFnegfDescr "-u"
double TFnegf(double u, double t);
double TFnegf(double u, double t)
  { return (-u); }
static TOFM TFnegfMap;

#define TFnhffDescr "-u/2"
double TFnhff(double u, double t);
double TFnhff(double u, double t)
  { return (-u/2); }
static TOFM TFnhffMap;

#define TFbumpDescr "exp(-t^2)"
double TFbump(double u, double t);
double TFbump(double u, double t)
  { return exp(-t*t); }
static TOFM TFbumpMap;

#define TFexntDescr "exp(-t)"
double TFexnt(double u, double t);
double TFexnt(double u, double t)
  { return exp(-t); }
static TOFM TFexntMap;

#define TFcostDescr "cos(t)"
double TFcost(double u, double t);
double TFcost(double u, double t)
  { return cos(t); }
static TOFM TFcostMap;

/* TABLE INITIALIZATION */

static bool_t SPTimeOnlyFuncMap_initialized = FALSE;

void SPTimeOnlyFuncMap_Initialize(void);
void SPTimeOnlyFuncMap_Initialize(void)
  {
    TFzeroMap =  (TOFM){&TFzero, FALSE, TFzeroDescr};
    TFstepMap =  (TOFM){&TFstep, FALSE, TFstepDescr};
    TFselfMap =  (TOFM){&TFself, FALSE, TFselfDescr};
    TFnegfMap =  (TOFM){&TFnegf, FALSE, TFnegfDescr};
    TFnhffMap =  (TOFM){&TFnhff, FALSE, TFnhffDescr};
    TFbumpMap =  (TOFM){&TFbump, FALSE, TFbumpDescr};
    TFexntMap =  (TOFM){&TFexnt, FALSE, TFexntDescr};
    TFcostMap =  (TOFM){&TFcost, FALSE, TFcostDescr};
    SPTimeOnlyFuncMap_initialized = TRUE;
  }

SPTimeOnlyFuncMap SPTimeOnlyFuncMap_FromName(char *name)
  {
    if (! SPTimeOnlyFuncMap_initialized) { SPTimeOnlyFuncMap_Initialize(); }
    if (strcmp(name, "TFzero") == 0) { return TFzeroMap; }
    if (strcmp(name, "TFstep") == 0) { return TFstepMap; }
    if (strcmp(name, "TFself") == 0) { return TFselfMap; }
    if (strcmp(name, "TFnegf") == 0) { return TFnegfMap; }
    if (strcmp(name, "TFnhff") == 0) { return TFnhffMap; }
    if (strcmp(name, "TFbump") == 0) { return TFbumpMap; }
    if (strcmp(name, "TFexnt") == 0) { return TFexntMap; }
    if (strcmp(name, "TFcost") == 0) { return TFcostMap; }
    affirm(FALSE, "unknown TimeOnlyFuncMap name");
    return NoTOFMap;
  }

