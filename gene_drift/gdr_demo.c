/* See {gdr_prob.h}  */
/* Last edited on 2023-06-01 07:28:34 by stolfi */

#define gdr_demo_C_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <cmp.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <argparser.h>
#include <vec.h> 

#include <gdr_demo.h>

/* INTERNAL PROTOTYPES */
    
/* EXTERNAL FUNCTIONS */

gdr_demo_parms_t *gdr_demo_parms_new(void)
  { 
    gdr_demo_parms_t *dmp = (gdr_demo_parms_t *)notnull(malloc(sizeof(gdr_demo_parms_t)), "no mem");
    
    *dmp = (gdr_demo_parms_t)
      { 

        .cMax =        -1,      /* Max children count. */
        .cAlpha =     NAN,      /* Decay ratio for child count distribution. */

        .fMin =        -1,      /* Minimum age to have children. */
        .fMax =        -1,      /* Maximum age to have children. */
      };

    return dmp;
  }

void gdr_demo_parms_free(gdr_demo_parms_t *dmp)
  { 
    free(dmp);
  }

double_vec_t gdr_demo_compute_child_count_distr(gdr_demo_parms_t *dmp)
  {
    int32_t cMax = dmp->cMax;
    demand(cMax >= 2, "invalid {cMax}");
    double_vec_t cProb = double_vec_new(cMax+1);
    
    /* Compute raw exponential distrib with {cProb.e[1]=1}, mean children count {cAvg}: */
    double p = 1.0;
    for (uint32_t c = 1;  c <= cMax; c++)
      { cProb.e[c] = p;
        p = p * dmp->cAlpha;
      }
      
    double unit = pow(0.1, dmp->cPrec);
      
    /* Hope that it converges: */
    for (uint32_t iter = 0;  iter < 5; iter++)
      { /* Compute mean child count {sum_dmp}: */
        double sum_dmp = 1.0e-200;
        for (uint32_t c = 1;  c <= cMax; c++) { sum_dmp += c*cProb.e[c]; }

        /* Rescale {cProb.e[1..cMax]} so that the mean is 1: */
        for (uint32_t c = 1;  c <= cMax; c++) { cProb.e[c] /= sum_dmp; }
        
        /* Round to multiple of {unit}, avoiding zero: */
        for (uint32_t c = 1;  c <= cMax; c++) 
          { double p = floor(cProb.e[c]/unit + 0.5)*unit;
            if (p == 0) { p = unit; }
            cProb.e[c] = p;
          }
      }
    
    /* Set {cProb[0]} so that the sum is 1: */
    double sum_p = 0;
    for (uint32_t c = 1;  c <= cMax; c++) { sum_p += cProb.e[c]; }   
    assert((sum_p > 0) && (sum_p < 1.0));
    cProb.e[0] = 1.0 - sum_p;
    
    double_vec_trim(&cProb, cMax+1);
    
    return cProb;
  }

void gdr_demo_show_parms(char *title, int32_t s, gdr_demo_parms_t *dmp)
  { fprintf(stderr, "%s", title);
    if (s >= 0) { fprintf(stderr, " -- sex %d", s); }
    fprintf(stderr, ":\n");

    fprintf(stderr, "  cMax = %d\n",     dmp->cMax);  
    fprintf(stderr, "  cAlpha = %.8f\n", dmp->cAlpha); 
    fprintf(stderr, "  fMin = %d\n",     dmp->fMin);   
    fprintf(stderr, "  fMax = %d\n",     dmp->fMax);   
  }
    
void gdr_demo_throw_children
  ( double_vec_t *cProb,
    int32_t fMin,
    int32_t fMax,
    int32_t *cNum_P, 
    int32_t cAge[]
  )
  {
    int32_t cMax = cProb->ne - 1;
    
    /* Choose the number of children according to {cProb}: */
    int32_t cNum = 0; /* Default if toss fails. */
    double coin = drandom();
    for (uint32_t c = 0;  c <= cMax; c++)
      { if (coin  < cProb->e[c]) { cNum = c; break; }
        coin = coin - cProb->e[c];
      }
    
    /* Choose the ages of the children: */
    int32_t aNum = fMax - fMin + 1; /* Num of ages to choose from. */
    if (cNum < aNum)
      { /* Choose distinct ages at random: */
        int32_t aDom[aNum];
        for (uint32_t ka = 0;  ka < aNum; ka++) { aDom[ka] = fMin + ka; }
        for (uint32_t kc = 0;  kc < cNum; kc++)
          { int32_t ka = int32_abrandom(kc, aNum-1);
            cAge[kc] = aDom[ka];
            if (ka != kc) { aDom[ka] = aDom[kc]; }
          }
      }
    else if (cNum > aNum)
      { /* Choose ages at random, don't care if distinct: */
        for (uint32_t kc = 0;  kc < cNum; kc++)
          { cAge[kc] = int32_abrandom(fMin, fMax); }
      }
    else /* {cNum == aNum} */
      { /* Choose one child per year: */
        for (uint32_t kc = 0;  kc < cNum; kc++)
          { cAge[kc] = fMin + kc; }
      }
    
    (*cNum_P) = cNum;
    
  }

void gdr_demo_write_child_distr_tex_table
  ( char *outPrefix, 
    char *tag,
    int32_t s,
    double_vec_t *cProb,
    int32_t cProbPrec,
    double_vec_t *cFreq,
    int32_t cFreqPrec
  )
  {
    char *fname = jsprintf("%s-%s-%d-probs.tex", outPrefix, tag, s);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    
    auto void wrpr(double_vec_t *pr, int32_t c, int32_t d);
      /* Prints " & {pr.e[c]}" to {wr} with {d} fraction digits,
        but writes "~" if it would print as zero. */
      
    /* Find {cMax}, the max {c} in all tables: */
    int32_t cMax = (cProb->ne > cFreq->ne ? cProb->ne : cFreq->ne) - 1;
    
    /* Compute the reverse cumulative distributions child counts: */
    double_vec_t cProbCum = gdr_demo_reverse_cumul_distr(cProb);
    double_vec_t cFreqCum = gdr_demo_reverse_cumul_distr(cFreq);

    fprintf(wr, "\\begin{tabular}{|r||r|r||r|r|}\n");
    fprintf(wr, "  \\hline\n");
    fprintf(wr, "  \\vstr ~ & \\multicolumn{4}{|c|}{%s} \\\\\n", (s == 0 ? "male" : "female"));
    fprintf(wr, "  \\cline{1-1} \\cline{2-3} \\cline{4-5}\n");
    fprintf(wr, "  \\vstr $c$ & $\\prb_%d(c)$ & $\\frq_%d(c)$ & $\\Prb_%d(c)$ & $\\Frq_%d(c)$ \\\\\n", s, s, s, s);
    fprintf(wr, "  \\cline{1-1} \\cline{2-3} \\cline{4-5}\n");
    for (uint32_t c = 0;  c <= cMax; c++) 
      { fprintf(wr,  "  %3d", c);
        wrpr(cProb, c, cProbPrec);
        wrpr(cFreq, c, cFreqPrec);
        wrpr(&(cProbCum), c, cProbPrec);
        wrpr(&(cFreqCum), c, cFreqPrec);
        fprintf(wr,  " \\\\\n");        
      }
    fprintf(wr, "  \\hline\n");
    fprintf(wr, "\\end{tabular}\n");

    fclose(wr);
    return;
    
    void wrpr(double_vec_t *cDistr, int32_t c, int32_t prec)
      { fprintf(wr, " &");
        double pc = (c < (int32_t)cDistr->ne ? cDistr->e[c] : 0.0);
        /* Format {pc} with {prec} decimal digits: */
        char *xpr = NULL;
        char *xpr = jsprintf("%*.*f", prec+2, prec, pc);

        /* Check if it prints as zero, if so print "~" instead: */
        char *q = xpr;
        while ((*q == ' ') || (*q == '0') || (*q == '.')) { q++; }
        fprintf(wr, " %*s", prec+2, (*q == 0 ? "~" : xpr));
        
        free(xpr);
      }
  }

double_vec_t gdr_demo_reverse_cumul_distr (double_vec_t *iProb)
  { int32_t iMax = iProb->ne - 1;
    /* Compute the complemented probabilities: */
    double_vec_t iProbCum =  double_vec_new(iMax + 1); /* Complemented cumulative probabilities */
    double pge = 0.0;
    for (int32_t i = iMax; i >= 0; i--)
      { pge = pge + iProb->e[i];
        iProbCum.e[i] = fmax(0.0, fmin(1.0, pge));
      }
    assert(fabs(pge) - 1.0 < 1.0e-10); /* allow for roundoff. */
    double_vec_trim(&iProbCum, iMax + 1);
    return iProbCum;
  }

void gdr_demo_exponential_distr(int32_t iMin, int32_t iMax, double alpha, double_vec_t *iProb)
  {
    demand((0 <= iMin) && (iMin <= iMax), "invalid {iMin,iMax}");
    double_vec_expand(iProb, iMax+1);
    for (uint32_t i = 0;  i < iMin; i++)
      { iProb->e[i] = 0.0; }
    double p = 1.0;
    for (int32_t i = iMin; i <= iMax; i++)
      { iProb->e[i] = p;
        p = alpha * p;
      }
    double_vec_trim(iProb, iMax+1);
  }

void gdr_demo_distr_normalize_sum(double val, double_vec_t *iProb)
  {
    demand(val > 0.0, "invalid target sum {val}");
    double sum = 0.0;
    for (uint32_t i = 0;  i < (int32_t)iProb->ne; i++) { sum += iProb->e[i]; }
    demand(sum > 0.0, "distribution is all zeros");
    double scale = val/sum;
    for (uint32_t i = 0;  i < (int32_t)iProb->ne; i++) { iProb->e[i] *= scale; }
  }

void gdr_demo_distr_round_off(double unit, double_vec_t *iProb)
  {
    bool_t debug = FALSE;
    
    for (uint32_t i = 0;  i < (int32_t)iProb->ne; i++)
      { double p = iProb->e[i];
        double p_round = floor(p/unit + 0.5)*unit;
        if (debug) { fprintf(stderr, "    raw %3d = %18.16f rounded %18.16f\n", i, p, p_round); } 
        iProb->e[i] = p_round;
      }
    /* Trim zero tail: */
    int32_t n = iProb->ne;
    while ((n > 0) && (iProb->e[n-1] == 0.0)) { n--; }
    double_vec_trim(iProb, n);
  }
    
void gdr_demo_distr_check(double_vec_t *cProb, int32_t cMax, double avg, int32_t prec)
  {
    demand(cProb->ne == cMax + 1, "inconsistent {cProb.ne,cMax}");
    double sum_p = 0;
    double sum_pc = 0;
    for (uint32_t c = 0;  c <= cMax; c++)
      { double p = cProb->e[c];
        demand((p >= 0.0) && (p <= 1.0), "invalid {cProb[c]}");
        sum_p += p;
        sum_pc += p*c;
      }
    demand(fabs(sum_p - 1.0) <= 2*pow(0.1, prec), "sum of {cProb} not 1");
    if (! isnan(avg))
      { demand(fabs(sum_pc - avg) <= cMax*pow(0.1, prec), "mean of {cProb} not 1"); }
  }

void gdr_demo_increment_count(int64_vec_t *rCount, int32_t *rMax_P, int32_t r, int32_t amt)
  {
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "          > %s\n", __FUNCTION__); }

    demand(r >= 0, "invalid index {r}");
    demand(amt >= 0, "invalid increment {amt}");
    if (amt > 0)
      { int32_t rMax = (*rMax_P);
        demand(rMax <= (int32_t)rCount->ne - 1, "inconsistent {rMax,ne}");
        int64_vec_expand(rCount, r);
        if (debug && (r > rMax)) { fprintf(stderr, "            expand from %d to %d\n", rMax, r); } 
        while (rMax < r) { rMax++; rCount->e[rMax] = 0; }
        rCount->e[r] += amt;
        if (debug) { fprintf(stderr, "            rCount[%d] = %ld\n", r, rCount->e[r]); } 
        (*rMax_P) = rMax;
      }
    if (debug) { fprintf(stderr, "          < %s\n", __FUNCTION__); }
  }

double_vec_t gdr_demo_freqs_from_counts(int32_t nr, int64_t rCount[])
  { bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "  > %s\n", __FUNCTION__); }

    /* Exclude the zero tail, if any: */
    while ((nr > 0) & (rCount[nr-1] == 0)) { nr--; }
    demand(nr > 0, "counts are all zero");
    
    /* Compute sum of entries: */
    int64_t sum = 0;
    for (uint32_t r = 0;  r < nr; r++) 
      { if (debug) { fprintf(stderr, "       r = %3d = %12ld\n", r, rCount[r]); }
        sum += rCount[r];
      }
    assert(sum > 0);
    
    /* Convert counts to frequencies: */
    double_vec_t rFreq = double_vec_new(nr);
    for (uint32_t r = 0;  r < nr; r++) 
      { rFreq.e[r] = ((double)(rCount[r]))/((double)sum); }
    double_vec_trim(&(rFreq), nr); /* Paranoia. */

    return rFreq;
  }

void gdr_demo_get_interesting_features(double_vec_t *cProb, int32_t *cBig_P, int32_t *cBigNum_P)
  { 
    /* Coose a remarkable number of sons: {cBigM} */
    int32_t cMax = cProb->ne - 1;
    int32_t cBig = cMax - 3;
    demand(cBig >= 2, "CD is too narrow (1)");
    /* Compute the probability of {cBig} sons or more: */
    double cBigProb = 0;
    for (int32_t c = cBig; c <= cMax; c++) 
      { cBigProb += cProb->e[c]; }
    int32_t cBigNum = (int32_t)floor(cBigProb*1000 + 1.0e-10);
    affirm(cBigNum >= 1, "CD is too narrow (2)");
    
    (*cBig_P) = cBig;
    (*cBigNum_P) = cBigNum;
  }

void gdr_demo_show_distr(char *title, int32_t s, double_vec_t *iProb, double_vec_t *iFreq)
  { fprintf(stderr, "%s", title);
    if (iFreq == NULL) 
      { fprintf(stderr, " distribution"); }
    else
      { fprintf(stderr, " probabilities and frequencies"); }
    if (s >= 0) { fprintf(stderr, " -- sex %d", s); }
    fprintf(stderr, ":\n");

    int32_t iMaxP = iProb->ne-1;
    demand(iProb->e[iMaxP] > 0.0, "trailing zero in {iProb}");
    int32_t iMaxF = (iFreq == NULL ? -1 : iFreq->ne-1);
    if (iFreq != NULL) { demand(iFreq->e[iMaxF] > 0.0, "trailing zero in {iFreq}"); }
    int32_t iMax = (int32_t)imax(iMaxP, iMaxF);

    auto void wrt(double_vec_t *vec, int32_t i);
      /* Prints {vec->e[i]} if {i < vec->ne}, else prints "". */
    
    auto void wrt_sum(double_vec_t *vec);
      /* Prints the sum of {vec->e[i]}. */

    auto void wrt_ave(double_vec_t *vec);
      /* Prints the mean of {vec->e[i]}. */
    
    fprintf(stderr, "  %3s %10s", "c/d", "prob");
    if (iFreq != NULL) { fprintf(stderr, " %10s", "freq"); }
    fprintf(stderr, "\n");

    fprintf(stderr, "  %3s %10s", "---", "--------");
    if (iFreq != NULL) { fprintf(stderr, " %10s", "--------"); }
    fprintf(stderr, "\n");

    for (uint32_t i = 0;  i <= iMax; i++)
      { fprintf(stderr, "  %3d", i);
        wrt(iProb, i);
        wrt(iFreq, i);
        fprintf(stderr, "\n");
      }

    fprintf(stderr, "\n");
   
    fprintf(stderr, "  %3s", "sum");
    wrt_sum(iProb);
    wrt_sum(iFreq);
    fprintf(stderr, "\n");
   
    fprintf(stderr, "  %3s", "ave");
    wrt_ave(iProb);
    wrt_ave(iFreq);
    fprintf(stderr, "\n");
      
    return;
    
    void wrt(double_vec_t *vec, int32_t i)
      { if (vec == NULL) { return; }
        if (i < (int32_t)vec->ne)
          { fprintf(stderr, " %10.6f", vec->e[i]); }
        else
          { fprintf(stderr, " %10s", ""); }
      }
      
    void wrt_sum(double_vec_t *vec)
      { if (vec == NULL) { return; }
        double sum = 0.0; 
        for (uint32_t i = 0;  i < (int32_t)vec->ne; i++) { sum += vec->e[i]; }
        fprintf(stderr, " %10.6f", sum);
      }
      
    void wrt_ave(double_vec_t *vec)
      { if (vec == NULL) { return; }
        double ave = 0.0; 
        for (uint32_t i = 0;  i < (int32_t)vec->ne; i++) { ave += i * vec->e[i]; }
        fprintf(stderr, " %10.6f", ave);
      }
  }

double gdr_demo_find_quantile(int32_t n, double x[], double frac)
  {
    demand(n > 0, "invalid {n}");
    demand((frac >= 0.0) && (frac <= 1.0), "invalid {frac}");
    double z = frac*n - 0.05;
    if (z <= 0)
      { return (double)x[0]; }
    else if (z >= n-1)
      { return (double)x[n-1]; }
    else
      { assert(n >= 2);
        int32_t k = (int32_t)floor(z);
        /* Just in case there is some funny rounding: */
        if (k >= n-1) { k = n-2; }
        assert((k >= 0) && (k < n));
        double r = z - k;
        double s = 1.0 - r;
        double v = s*x[k] + r*x[k+1];
        return v;
      }
  }

void gdr_demo_sort_doubles(int32_t n, double x[])
  { 
    auto int cmp(const void *a, const void *b);
    
    qsort(x, n, sizeof(double), &cmp);
    /* Paranoia: */
    for (uint32_t r = 1;  r < n; r++) { assert(x[r-1] <= x[r]); }
    return;
    
    int cmp(const void *a, const void *b)
      { double *ai = (double *)a;
        double *bi = (double *)b;
        return cmp_double(ai, bi);
      }
  }

gdr_demo_parms_t *gdr_demo_options_parse(argparser_t *pp)
  { gdr_demo_parms_t *dmp = gdr_demo_parms_new();
    
    dmp->cMax =   (int32_t)argparser_get_next_int(pp, 2, gdr_demo_child_count_MAX); 
    dmp->cAlpha = argparser_get_next_double(pp, 0.05, 0.99);
    dmp->cPrec =  (int32_t)argparser_get_next_int(pp, 1, 6);
    dmp->fMin =   (int32_t)argparser_get_next_int(pp, 1, gdr_demo_age_MAX); 
    dmp->fMax =   (int32_t)argparser_get_next_int(pp, dmp->fMin, gdr_demo_age_MAX);  
    
    return dmp;
  }
