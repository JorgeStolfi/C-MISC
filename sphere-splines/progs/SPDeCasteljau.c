/* See SPDeCasteljau.h */
/* Last edited on 2005-10-27 17:07:11 by stolfi */

#include <SPDeCasteljau.h>
#include <SPBasic.h>

#include <r3.h>
#include <vec.h>
#include <affirm.h>
#include <nat.h>

BezLabel SPDeCasteljau_IndexToBezLabel(nat_t t, nat_t deg)
  { nat_t r = 0, s = 0, k = t;
    while (k > r) { r = r + 1; s = s + r; k = k - r; }
    return (BezLabel){{deg - r, s + r - t, t - s}};
  }

nat_t SPDeCasteljau_BezLabelToIndex(BezLabel e, nat_t deg)
  { /* Should be made faster... */
    nat_t r, s = 0, i = e.e[0], k = e.e[2];
    for (r = i; r < deg; r++) { s = s + (deg-r); }
    nat_t t = s + k;
    return t;
  }

void SPDeCasteljau_EvalGen(BezCoeff_vec_t coeff, int deg, r3_t *a, double *f)
  { int nc = (deg + 1)*(deg + 2)/2;
    affirm(coeff.ne == nc, "wrong num of coeffs");
    double a0 = a->c[0], a1 = a->c[1], a2 = a->c[2]; 
    if (deg == 0)
      { (*f) = coeff.e[0]; }
    else if (deg == 1)
      { (*f) = coeff.e[0]*a0 + coeff.e[1]*a1 + coeff.e[2]*a2; }
    else
      { double c[nc];
        /* {c[R]} is coefficient {R} of an intermediate spline. */
        /* Copy coeffs to work coeffs {c}: */
        { nat_t r; 
          for (r = 0; r < nc; r++)
           { c[r] = coeff.e[r]; }
        }
        /* Apply the DeCasteljau recursion: */
        while(deg > 0)
          { int degn = deg-1;
            int i, j, k;
            nat_t r0 = 0;
            for (i = degn ; i >= 0; i--)
              { for (j = degn - i; j >= 0; j--)
                  { k = degn - i - j;
                    nat_t r1 = r0 + (degn - i + 1);
                    nat_t r2 = r1 + 1;
                    c[r0] = a0*c[r0] + a1*c[r1] + a2*c[r2];
                    r0++;
                  }
              }
            deg = degn;
          }
        (*f) = c[0];
      }
  }

void SPDeCasteljau_GradGen(BezCoeff_vec_t coeff, int deg, r3_t *a, double *f, r3_t *df)
  { int nc = (deg + 1)*(deg + 2)/2;
    affirm(coeff.ne == nc, "wrong num of coeffs");
    double a0 = a->c[0], a1 = a->c[1], a2 = a->c[2]; 
    if (deg == 0)
      { (*f) = coeff.e[0];
        (*df) = (r3_t){{ 0.0, 0.0, 0.0 }};
      }
    else if (deg == 1)
      { (*f) = coeff.e[0]*a0 + coeff.e[1]*a1 + coeff.e[2]*a2;
        (*df) = (r3_t){{ coeff.e[0], coeff.e[1], coeff.e[2] }};
      }
    else
      { double c[nc]; 
        /* {c[R]} is coefficient {R} of an intermediate spline. */
        double d0[nc], d1[nc], d2[nc];
        /* {dN[R]} is the derivative of {c[R]} with respect to {a[N]}. */
        /* Copy coeffs to work coeffs {c}, intialize derivatives: */
        { nat_t r; 
          for (r = 0; r < nc; r++) 
            { c[r] = coeff.e[r]; d0[r] = d1[r] = d2[r] = 0.0; }
        }
        /* Apply the DeCasteljau recursion: */
        while(deg > 0)
          { int degn = deg-1;
            int i, j, k;
            nat_t r0 = 0;
            for (i = degn ; i >= 0; i--)
              { for (j = degn - i; j >= 0; j--)
                  { k = degn - i - j;
                    nat_t r1 = r0 + (degn - i + 1);
                    nat_t r2 = r1 + 1;
                    d0[r0] = c[r0] + a0*d0[r0] + a1*d0[r1] + a2*d0[r2];
                    d1[r0] = c[r1] + a0*d1[r0] + a1*d1[r1] + a2*d1[r2];
                    d2[r0] = c[r2] + a0*d2[r0] + a1*d2[r1] + a2*d2[r2];
                    c[r0] = a0*c[r0] + a1*c[r1] + a2*c[r2];
                    r0++;
                  }
              }
            deg = degn;
          }
        (*f) = c[0];
        (*df) = (r3_t){{d0[0], d1[0], d2[0]}};
      }
  }

void SPDeCasteljau_HessGen(BezCoeff_vec_t coeff, int deg, r3_t *a, double *f, r3_t *df, r6_t *ddf)
  { int nc = (deg + 1)*(deg + 2)/2;
    affirm(coeff.ne == nc, "wrong num of coeffs");
    double a0 = a->c[0], a1 = a->c[1], a2 = a->c[2]; 
    if (deg == 0)
      { (*f) = coeff.e[0];
        (*df) = (r3_t){{ 0.0, 0.0, 0.0 }};
        (*ddf) = (r6_t){{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
      }
    else if (deg == 1)
      { (*f) = coeff.e[0]*a0 + coeff.e[1]*a1 + coeff.e[2]*a2;
        (*df) = (r3_t){{ coeff.e[0], coeff.e[1], coeff.e[2] }};
        (*ddf) = (r6_t){{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }};
      }
    else
      { double c[nc]; 
        /* {c[R]} is coefficient {R} of an intermediate spline. */
        double d0[nc], d1[nc], d2[nc];
        /* {dN[R]} is the derivative of {c[R]} with respect to {a[N]}. */
        double d00[nc], d01[nc], d11[nc], d02[nc], d12[nc], d22[nc];
        /* {dMN[R]} is the derivative of {c[R]} with respect to {a[M]}  and {a[N]}. */
        /* Copy coeffs to work coeffs {c}, intialize derivatives: */
        { nat_t r; 
          for (r = 0; r < nc; r++) 
            { c[r] = coeff.e[r];
              d0[r] = d1[r] = d2[r] = 0.0;
              d00[r] = d01[r] = d11[r] = d02[r] = d12[r] = d22[r] = 0.0;
            }
        }
        /* Apply the DeCasteljau recursion: */
        while(deg > 0)
          { int degn = deg-1;
            int i, j, k;
            nat_t r0 = 0;
            for (i = degn ; i >= 0; i--)
              { for (j = degn - i; j >= 0; j--)
                  { k = degn - i - j;
                    nat_t r1 = r0 + (degn - i + 1);
                    nat_t r2 = r1 + 1;

                    d00[r0] = d0[r0] + d0[r0] + a0*d00[r0] + a1*d00[r1] + a2*d00[r2];
                    d01[r0] = d1[r0] + d0[r1] + a0*d01[r0] + a1*d01[r1] + a2*d01[r2];
                    d11[r0] = d1[r1] + d1[r1] + a0*d11[r0] + a1*d11[r1] + a2*d11[r2];
                    d02[r0] = d2[r0] + d0[r2] + a0*d02[r0] + a1*d02[r1] + a2*d02[r2];
                    d12[r0] = d2[r1] + d1[r2] + a0*d12[r0] + a1*d12[r1] + a2*d12[r2];
                    d22[r0] = d2[r2] + d2[r2] + a0*d22[r0] + a1*d22[r1] + a2*d22[r2];

                    d0[r0] = c[r0] + a0*d0[r0] + a1*d0[r1] + a2*d0[r2];
                    d1[r0] = c[r1] + a0*d1[r0] + a1*d1[r1] + a2*d1[r2];
                    d2[r0] = c[r2] + a0*d2[r0] + a1*d2[r1] + a2*d2[r2];

                    c[r0] = a0*c[r0] + a1*c[r1] + a2*c[r2];
                    r0++;
                  }
              }
            deg = degn;
          }
        (*f) = c[0];
        (*df) = (r3_t){{d0[0], d1[0], d2[0]}};
        (*ddf) = (r6_t){{d00[0], d01[0], d11[0], d02[0], d12[0], d22[0]}};
      }
  }

/* SPECIALIZED VERSIONS */

void SPDeCasteljau_Eval0(BezCoeff_vec_t coeff, r3_t *a, double *f)
  {
    affirm(coeff.ne == 1, "wrong num of coeffs");
    
    double c000 = coeff.e[0];
    
    (*f) = c000;
  }

void SPDeCasteljau_Eval1(BezCoeff_vec_t coeff, r3_t *a, double *f)
  {
    affirm(coeff.ne == 3, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c100 = coeff.e[0];
    double c010 = coeff.e[1];
    double c001 = coeff.e[2];
          
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;

    (*f) = c000;
  }

void SPDeCasteljau_Eval2(BezCoeff_vec_t coeff, r3_t *a, double *f)
  {
    affirm(coeff.ne == 6, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c200 = coeff.e[0];
    double c110 = coeff.e[1];
    double c101 = coeff.e[2];
    double c020 = coeff.e[3];
    double c011 = coeff.e[4];
    double c002 = coeff.e[5];
          
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;

    (*f) = c000;
  }

void SPDeCasteljau_Eval3(BezCoeff_vec_t coeff, r3_t *a, double *f)
  {
    affirm(coeff.ne = 10, "wrong num of coefs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c300 = coeff.e[0];
    double c210 = coeff.e[1];
    double c201 = coeff.e[2];
    double c120 = coeff.e[3];
    double c111 = coeff.e[4];
    double c102 = coeff.e[5];
    double c030 = coeff.e[6];
    double c021 = coeff.e[7];
    double c012 = coeff.e[8];
    double c003 = coeff.e[9];
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;

    (*f) = c000;
  }
  
void SPDeCasteljau_Eval4(BezCoeff_vec_t coeff, r3_t *a, double *f) 
  {
    affirm(coeff.ne == 15, "wrong num of coeffs");

    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c400 = coeff.e[ 0];
    double c310 = coeff.e[ 1];
    double c301 = coeff.e[ 2];
    double c220 = coeff.e[ 3];
    double c211 = coeff.e[ 4];
    double c202 = coeff.e[ 5];
    double c130 = coeff.e[ 6];
    double c121 = coeff.e[ 7];
    double c112 = coeff.e[ 8];
    double c103 = coeff.e[ 9];
    double c040 = coeff.e[10];
    double c031 = coeff.e[11];
    double c022 = coeff.e[12];
    double c013 = coeff.e[13];
    double c004 = coeff.e[14];
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;

    (*f) = c000;
  }
  
void SPDeCasteljau_Eval5(BezCoeff_vec_t coeff, r3_t *a, double *f)
  {
    affirm(coeff.ne == 21, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c500 = coeff.e[ 0];
    double c410 = coeff.e[ 1];
    double c401 = coeff.e[ 2];
    double c320 = coeff.e[ 3];
    double c311 = coeff.e[ 4];
    double c302 = coeff.e[ 5];
    double c230 = coeff.e[ 6];
    double c221 = coeff.e[ 7];
    double c212 = coeff.e[ 8];
    double c203 = coeff.e[ 9];
    double c140 = coeff.e[10];
    double c131 = coeff.e[11];
    double c122 = coeff.e[12];
    double c113 = coeff.e[13];
    double c104 = coeff.e[14];
    double c050 = coeff.e[15];
    double c041 = coeff.e[16];
    double c032 = coeff.e[17];
    double c023 = coeff.e[18];
    double c014 = coeff.e[19];
    double c005 = coeff.e[20];
    
    double c400 = a0 * c500 + a1 * c410 + a2 * c401;
    double c310 = a0 * c410 + a1 * c320 + a2 * c311;
    double c301 = a0 * c401 + a1 * c311 + a2 * c302;
    double c220 = a0 * c320 + a1 * c230 + a2 * c221;
    double c211 = a0 * c311 + a1 * c221 + a2 * c212;
    double c202 = a0 * c302 + a1 * c212 + a2 * c203;
    double c112 = a0 * c212 + a1 * c122 + a2 * c113;
    double c121 = a0 * c221 + a1 * c131 + a2 * c122;
    double c103 = a0 * c203 + a1 * c113 + a2 * c104;
    double c130 = a0 * c230 + a1 * c140 + a2 * c131;
    double c022 = a0 * c122 + a1 * c032 + a2 * c023;
    double c013 = a0 * c113 + a1 * c023 + a2 * c014;
    double c031 = a0 * c131 + a1 * c041 + a2 * c032;
    double c040 = a0 * c140 + a1 * c050 + a2 * c041;
    double c004 = a0 * c104 + a1 * c014 + a2 * c005;
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;

    (*f) = c000;
  }

   
void SPDeCasteljau_Eval6(BezCoeff_vec_t coeff, r3_t *a, double *f)
  {
    affirm(coeff.ne == 28, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c600 = coeff.e[ 0];
    double c510 = coeff.e[ 1];
    double c501 = coeff.e[ 2];
    double c420 = coeff.e[ 3];
    double c411 = coeff.e[ 4];
    double c402 = coeff.e[ 5];
    double c330 = coeff.e[ 6];
    double c321 = coeff.e[ 7];
    double c312 = coeff.e[ 8];
    double c303 = coeff.e[ 9];
    double c240 = coeff.e[10];
    double c231 = coeff.e[11];
    double c222 = coeff.e[12];
    double c213 = coeff.e[13];
    double c204 = coeff.e[14];
    double c150 = coeff.e[15];
    double c141 = coeff.e[16];
    double c132 = coeff.e[17];
    double c123 = coeff.e[18];
    double c114 = coeff.e[19];
    double c105 = coeff.e[20];
    double c060 = coeff.e[21];
    double c051 = coeff.e[22];
    double c042 = coeff.e[23];
    double c033 = coeff.e[24];
    double c024 = coeff.e[25];
    double c015 = coeff.e[26];
    double c006 = coeff.e[27];
    
    double c500 = a0 * c600 + a1 * c510 + a2 * c501;
    double c410 = a0 * c510 + a1 * c420 + a2 * c411;
    double c401 = a0 * c501 + a1 * c411 + a2 * c402;
    double c320 = a0 * c420 + a1 * c330 + a2 * c321;
    double c311 = a0 * c411 + a1 * c321 + a2 * c312;
    double c302 = a0 * c402 + a1 * c312 + a2 * c303;
    double c230 = a0 * c330 + a1 * c240 + a2 * c231;
    double c221 = a0 * c321 + a1 * c231 + a2 * c222;
    double c212 = a0 * c312 + a1 * c222 + a2 * c213;
    double c203 = a0 * c303 + a1 * c213 + a2 * c204;
    double c140 = a0 * c240 + a1 * c150 + a2 * c141;
    double c131 = a0 * c231 + a1 * c141 + a2 * c132;
    double c122 = a0 * c222 + a1 * c132 + a2 * c123;
    double c113 = a0 * c213 + a1 * c123 + a2 * c114;
    double c104 = a0 * c204 + a1 * c114 + a2 * c105;
    double c050 = a0 * c150 + a1 * c060 + a2 * c051;
    double c041 = a0 * c141 + a1 * c051 + a2 * c042;
    double c032 = a0 * c132 + a1 * c042 + a2 * c033;
    double c023 = a0 * c123 + a1 * c033 + a2 * c024;
    double c014 = a0 * c114 + a1 * c024 + a2 * c015;
    double c005 = a0 * c105 + a1 * c015 + a2 * c006;
    
    
    double c400 = a0 * c500 + a1 * c410 + a2 * c401;
    double c310 = a0 * c410 + a1 * c320 + a2 * c311;
    double c301 = a0 * c401 + a1 * c311 + a2 * c302;
    double c220 = a0 * c320 + a1 * c230 + a2 * c221;
    double c211 = a0 * c311 + a1 * c221 + a2 * c212;
    double c202 = a0 * c302 + a1 * c212 + a2 * c203;
    double c112 = a0 * c212 + a1 * c122 + a2 * c113;
    double c121 = a0 * c221 + a1 * c131 + a2 * c122;
    double c103 = a0 * c203 + a1 * c113 + a2 * c104;
    double c130 = a0 * c230 + a1 * c140 + a2 * c131;
    double c022 = a0 * c122 + a1 * c032 + a2 * c023;
    double c013 = a0 * c113 + a1 * c023 + a2 * c014;
    double c031 = a0 * c131 + a1 * c041 + a2 * c032;
    double c040 = a0 * c140 + a1 * c050 + a2 * c041;
    double c004 = a0 * c104 + a1 * c014 + a2 * c005;
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;

    (*f) = c000;
  }

void SPDeCasteljau_Eval7(BezCoeff_vec_t coeff, r3_t *a, double *f)
  {
    affirm(coeff.ne == 36, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c700 = coeff.e[ 0];
    double c610 = coeff.e[ 1];
    double c601 = coeff.e[ 2];
    double c520 = coeff.e[ 3];
    double c511 = coeff.e[ 4];
    double c502 = coeff.e[ 5];
    double c430 = coeff.e[ 6];
    double c421 = coeff.e[ 7];
    double c412 = coeff.e[ 8];
    double c403 = coeff.e[ 9];
    double c340 = coeff.e[10];
    double c331 = coeff.e[11];
    double c322 = coeff.e[12];
    double c313 = coeff.e[13];
    double c304 = coeff.e[14];
    double c250 = coeff.e[15];
    double c241 = coeff.e[16];
    double c232 = coeff.e[17];
    double c223 = coeff.e[18];
    double c214 = coeff.e[19];
    double c205 = coeff.e[20];
    double c160 = coeff.e[21];
    double c151 = coeff.e[22];
    double c142 = coeff.e[23];
    double c133 = coeff.e[24];
    double c124 = coeff.e[25];
    double c115 = coeff.e[26];
    double c106 = coeff.e[27];
    double c070 = coeff.e[28];
    double c061 = coeff.e[29];
    double c052 = coeff.e[30];
    double c043 = coeff.e[31];
    double c034 = coeff.e[32];
    double c025 = coeff.e[33];
    double c016 = coeff.e[34];
    double c007 = coeff.e[35];
    
    double c600 = a0 * c700 + a1 * c610 + a2 * c601;
    double c510 = a0 * c610 + a1 * c520 + a2 * c511;
    double c501 = a0 * c601 + a1 * c511 + a2 * c502;
    double c420 = a0 * c520 + a1 * c430 + a2 * c421;
    double c411 = a0 * c511 + a1 * c421 + a2 * c412;
    double c402 = a0 * c502 + a1 * c412 + a2 * c403;
    double c330 = a0 * c430 + a1 * c340 + a2 * c331;
    double c321 = a0 * c421 + a1 * c331 + a2 * c322;
    double c312 = a0 * c412 + a1 * c322 + a2 * c313;
    double c303 = a0 * c403 + a1 * c313 + a2 * c304;
    double c240 = a0 * c340 + a1 * c250 + a2 * c241;
    double c231 = a0 * c331 + a1 * c241 + a2 * c232;
    double c222 = a0 * c322 + a1 * c232 + a2 * c223;
    double c213 = a0 * c313 + a1 * c223 + a2 * c214;
    double c204 = a0 * c304 + a1 * c214 + a2 * c205;
    double c150 = a0 * c250 + a1 * c160 + a2 * c151;
    double c141 = a0 * c241 + a1 * c151 + a2 * c142;
    double c132 = a0 * c232 + a1 * c142 + a2 * c133;
    double c123 = a0 * c223 + a1 * c133 + a2 * c124;
    double c114 = a0 * c214 + a1 * c124 + a2 * c115;
    double c105 = a0 * c205 + a1 * c115 + a2 * c106;
    double c060 = a0 * c160 + a1 * c070 + a2 * c061;
    double c051 = a0 * c151 + a1 * c061 + a2 * c052;
    double c042 = a0 * c142 + a1 * c052 + a2 * c043;
    double c033 = a0 * c133 + a1 * c043 + a2 * c034;
    double c024 = a0 * c124 + a1 * c034 + a2 * c025;
    double c015 = a0 * c115 + a1 * c025 + a2 * c016;
    double c006 = a0 * c106 + a1 * c016 + a2 * c007;

    double c500 = a0 * c600 + a1 * c510 + a2 * c501;
    double c410 = a0 * c510 + a1 * c420 + a2 * c411;
    double c401 = a0 * c501 + a1 * c411 + a2 * c402;
    double c320 = a0 * c420 + a1 * c330 + a2 * c321;
    double c311 = a0 * c411 + a1 * c321 + a2 * c312;
    double c302 = a0 * c402 + a1 * c312 + a2 * c303;
    double c230 = a0 * c330 + a1 * c240 + a2 * c231;
    double c221 = a0 * c321 + a1 * c231 + a2 * c222;
    double c212 = a0 * c312 + a1 * c222 + a2 * c213;
    double c203 = a0 * c303 + a1 * c213 + a2 * c204;
    double c140 = a0 * c240 + a1 * c150 + a2 * c141;
    double c131 = a0 * c231 + a1 * c141 + a2 * c132;
    double c122 = a0 * c222 + a1 * c132 + a2 * c123;
    double c113 = a0 * c213 + a1 * c123 + a2 * c114;
    double c104 = a0 * c204 + a1 * c114 + a2 * c105;
    double c050 = a0 * c150 + a1 * c060 + a2 * c051;
    double c041 = a0 * c141 + a1 * c051 + a2 * c042;
    double c032 = a0 * c132 + a1 * c042 + a2 * c033;
    double c023 = a0 * c123 + a1 * c033 + a2 * c024;
    double c014 = a0 * c114 + a1 * c024 + a2 * c015;
    double c005 = a0 * c105 + a1 * c015 + a2 * c006;
    
    double c400 = a0 * c500 + a1 * c410 + a2 * c401;
    double c310 = a0 * c410 + a1 * c320 + a2 * c311;
    double c301 = a0 * c401 + a1 * c311 + a2 * c302;
    double c220 = a0 * c320 + a1 * c230 + a2 * c221;
    double c211 = a0 * c311 + a1 * c221 + a2 * c212;
    double c202 = a0 * c302 + a1 * c212 + a2 * c203;
    double c130 = a0 * c230 + a1 * c140 + a2 * c131;
    double c121 = a0 * c221 + a1 * c131 + a2 * c122;
    double c112 = a0 * c212 + a1 * c122 + a2 * c113;
    double c103 = a0 * c203 + a1 * c113 + a2 * c104;
    double c040 = a0 * c140 + a1 * c050 + a2 * c041;
    double c031 = a0 * c131 + a1 * c041 + a2 * c032;
    double c022 = a0 * c122 + a1 * c032 + a2 * c023;
    double c013 = a0 * c113 + a1 * c023 + a2 * c014;
    double c004 = a0 * c104 + a1 * c014 + a2 * c005;
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;

    (*f) = c000;
  }

void SPDeCasteljau_Grad0(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df)
  {
    affirm(coeff.ne == 1, "wrong num of coeffs");
    
    double c000 = coeff.e[0];

    (*f) = c000;
    (*df) = (r3_t){{0.0, 0.0, 0.0}};
  }

void SPDeCasteljau_Grad1(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df)
  {
    affirm(coeff.ne == 3, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c100 = coeff.e[0];
    double c010 = coeff.e[1];
    double c001 = coeff.e[2];
          
    /* dXYZdK is the derivative of cXYZ with respect to aK. */
    
    /* Layer 0 */
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;
    double d000d0 = c100;
    double d000d1 = c010;
    double d000d2 = c001;

    (*f) = c000;
    (*df) = (r3_t){{d000d0, d000d1, d000d2}};
    
  }

void SPDeCasteljau_Grad2(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df)
  {
    affirm(coeff.ne == 6, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c200 = coeff.e[0];
    double c110 = coeff.e[1];
    double c101 = coeff.e[2];
    double c020 = coeff.e[3];
    double c011 = coeff.e[4];
    double c002 = coeff.e[5];
          
    /* dXYZdK is the derivative of cXYZ with respect to aK. */
    
    /* Layer 1 */
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double d100d0 = c200;
    double d100d1 = c110;
    double d100d2 = c101;
         
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double d010d0 = c110;
    double d010d1 = c020;
    double d010d2 = c011;
    
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    double d001d0 = c101;
    double d001d1 = c011;
    double d001d2 = c002;
    
    /* Layer 0 */
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;
    double d000d0 = c100 + a0*d100d0 + a1*d010d0 + a2*d001d0;
    double d000d1 = c010 + a0*d100d1 + a1*d010d1 + a2*d001d1;
    double d000d2 = c001 + a0*d100d2 + a1*d010d2 + a2*d001d2;

    (*f) = c000;
    (*df) = (r3_t){{d000d0, d000d1, d000d2}};
  }

void SPDeCasteljau_Grad3(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df)
  {
    affirm(coeff.ne == 10, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c300 = coeff.e[0];
    double c210 = coeff.e[1];
    double c201 = coeff.e[2];
    double c120 = coeff.e[3];
    double c111 = coeff.e[4];
    double c102 = coeff.e[5];
    double c030 = coeff.e[6];
    double c021 = coeff.e[7];
    double c012 = coeff.e[8];
    double c003 = coeff.e[9];
    
    /* dXYZdK is the derivative of cXYZ with respect to aK. */
    
    /* Layer 2 */  
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double d200d0 = c300;
    double d200d1 = c210;
    double d200d2 = c201;
         
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double d110d0 = c210;
    double d110d1 = c120;
    double d110d2 = c111;
    
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double d101d0 = c201;
    double d101d1 = c111;
    double d101d2 = c102;
    
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double d020d0 = c120;
    double d020d1 = c030;
    double d020d2 = c021;
    
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double d011d0 = c111;
    double d011d1 = c021;
    double d011d2 = c012;
    
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    double d002d0 = c102;
    double d002d1 = c012;
    double d002d2 = c003;
    
    /* Layer 1 */
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double d100d0 = c200 + a0*d200d0 + a1*d110d0 + a2*d101d0;
    double d100d1 = c110 + a0*d200d1 + a1*d110d1 + a2*d101d1;
    double d100d2 = c101 + a0*d200d2 + a1*d110d2 + a2*d101d2;
          
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double d010d0 = c110 + a0*d110d0 + a1*d020d0 + a2*d011d0;
    double d010d1 = c020 + a0*d110d1 + a1*d020d1 + a2*d011d1;
    double d010d2 = c011 + a0*d110d2 + a1*d020d2 + a2*d011d2;
    
    
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    double d001d0 = c101 + a0*d101d0 + a1*d011d0 + a2*d002d0;
    double d001d1 = c011 + a0*d101d1 + a1*d011d1 + a2*d002d1;
    double d001d2 = c002 + a0*d101d2 + a1*d011d2 + a2*d002d2;

    /* Layer 0 */
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001; 
    double d000d0 = c100 + a0*d100d0 + a1*d010d0 + a2*d001d0;
    double d000d1 = c010 + a0*d100d1 + a1*d010d1 + a2*d001d1;
    double d000d2 = c001 + a0*d100d2 + a1*d010d2 + a2*d001d2;

    (*f) = c000;
    (*df) = (r3_t){{d000d0, d000d1, d000d2}};
  }

void SPDeCasteljau_Grad4(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df) 
  {
    affirm(coeff.ne == 15, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c400 = coeff.e[ 0];
    double c310 = coeff.e[ 1];
    double c301 = coeff.e[ 2];
    double c220 = coeff.e[ 3];
    double c211 = coeff.e[ 4];
    double c202 = coeff.e[ 5];
    double c130 = coeff.e[ 6];
    double c121 = coeff.e[ 7];
    double c112 = coeff.e[ 8];
    double c103 = coeff.e[ 9];
    double c040 = coeff.e[10];
    double c031 = coeff.e[11];
    double c022 = coeff.e[12];
    double c013 = coeff.e[13];
    double c004 = coeff.e[14];
    
    /* dXYZdK is the derivative of cXYZ with respect to aK. */
    
    /* Layer 3 */
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double d300d0 = c400;
    double d300d1 = c310;
    double d300d2 = c301;
    
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double d210d0 = c310;
    double d210d1 = c220;
    double d210d2 = c211;
    
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double d201d0 = c301;
    double d201d1 = c211;
    double d201d2 = c202;
    
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double d120d0 = c220;
    double d120d1 = c130;
    double d120d2 = c121;
    
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double d111d0 = c211;
    double d111d1 = c121;
    double d111d2 = c112;
    
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double d102d0 = c202;
    double d102d1 = c112;
    double d102d2 = c103;
    
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double d030d0 = c130;
    double d030d1 = c040;
    double d030d2 = c031;
    
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double d021d0 = c121;
    double d021d1 = c031;
    double d021d2 = c022;
    
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double d012d0 = c112;
    double d012d1 = c022;
    double d012d2 = c013;
    
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    double d003d0 = c103;
    double d003d1 = c013;
    double d003d2 = c004;
         
    /* Layer 2 */
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double d200d0 = c300 + a0*d300d0 + a1*d210d0 + a2*d201d0;
    double d200d1 = c210 + a0*d300d1 + a1*d210d1 + a2*d201d1;
    double d200d2 = c201 + a0*d300d2 + a1*d210d2 + a2*d201d2;
    
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double d110d0 = c210 + a0*d210d0 + a1*d120d0 + a2*d111d0;
    double d110d1 = c120 + a0*d210d1 + a1*d120d1 + a2*d111d1;
    double d110d2 = c111 + a0*d210d2 + a1*d120d2 + a2*d111d2;
    
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double d101d0 = c201 + a0*d201d0 + a1*d111d0 + a2*d102d0;
    double d101d1 = c111 + a0*d201d1 + a1*d111d1 + a2*d102d1;
    double d101d2 = c102 + a0*d201d2 + a1*d111d2 + a2*d102d2;
    
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double d020d0 = c120 + a0*d120d0 + a1*d030d0 + a2*d021d0;
    double d020d1 = c030 + a0*d120d1 + a1*d030d1 + a2*d021d1;
    double d020d2 = c021 + a0*d120d2 + a1*d030d2 + a2*d021d2;
    
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double d011d0 = c111 + a0*d111d0 + a1*d021d0 + a2*d012d0;
    double d011d1 = c021 + a0*d111d1 + a1*d021d1 + a2*d012d1;
    double d011d2 = c012 + a0*d111d2 + a1*d021d2 + a2*d012d2;
    
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    double d002d0 = c102 + a0*d102d0 + a1*d012d0 + a2*d003d0;
    double d002d1 = c012 + a0*d102d1 + a1*d012d1 + a2*d003d1;
    double d002d2 = c003 + a0*d102d2 + a1*d012d2 + a2*d003d2;
    
    /* Layer 1 */
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double d100d0 = c200 + a0*d200d0 + a1*d110d0 + a2*d101d0;
    double d100d1 = c110 + a0*d200d1 + a1*d110d1 + a2*d101d1;
    double d100d2 = c101 + a0*d200d2 + a1*d110d2 + a2*d101d2;
          
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double d010d0 = c110 + a0*d110d0 + a1*d020d0 + a2*d011d0;
    double d010d1 = c020 + a0*d110d1 + a1*d020d1 + a2*d011d1;
    double d010d2 = c011 + a0*d110d2 + a1*d020d2 + a2*d011d2;
    
    
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    double d001d0 = c101 + a0*d101d0 + a1*d011d0 + a2*d002d0;
    double d001d1 = c011 + a0*d101d1 + a1*d011d1 + a2*d002d1;
    double d001d2 = c002 + a0*d101d2 + a1*d011d2 + a2*d002d2;
    
    /* Layer 0 */
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;
    double d000d0 = c100 + a0*d100d0 + a1*d010d0 + a2*d001d0;
    double d000d1 = c010 + a0*d100d1 + a1*d010d1 + a2*d001d1;
    double d000d2 = c001 + a0*d100d2 + a1*d010d2 + a2*d001d2;

    (*f) = c000;
    (*df) = (r3_t){{d000d0, d000d1, d000d2}};
  }
  
void SPDeCasteljau_Grad5(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df)
   {
    affirm(coeff.ne == 21, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c500 = coeff.e[ 0];
    double c410 = coeff.e[ 1];
    double c401 = coeff.e[ 2];
    double c320 = coeff.e[ 3];
    double c311 = coeff.e[ 4];
    double c302 = coeff.e[ 5];
    double c230 = coeff.e[ 6];
    double c221 = coeff.e[ 7];
    double c212 = coeff.e[ 8];
    double c203 = coeff.e[ 9];
    double c140 = coeff.e[10];
    double c131 = coeff.e[11];
    double c122 = coeff.e[12];
    double c113 = coeff.e[13];
    double c104 = coeff.e[14];
    double c050 = coeff.e[15];
    double c041 = coeff.e[16];
    double c032 = coeff.e[17];
    double c023 = coeff.e[18];
    double c014 = coeff.e[19];
    double c005 = coeff.e[20];

    /* dXYZdK is the derivative of cXYZ with respect to aK. */
    
    /* Layer 4 */
    
    double c400 = a0 * c500 + a1 * c410 + a2 * c401;
    double d400d0 = c500;
    double d400d1 = c410;
    double d400d2 = c401;

    double c310 = a0 * c410 + a1 * c320 + a2 * c311;
    double d310d0 = c410;
    double d310d1 = c320;
    double d310d2 = c311;
    
    double c301 = a0 * c401 + a1 * c311 + a2 * c302;
    double d301d0 = c401;
    double d301d1 = c311;
    double d301d2 = c302;
    
    double c220 = a0 * c320 + a1 * c230 + a2 * c221;
    double d220d0 = c320;
    double d220d1 = c230;
    double d220d2 = c221;
    
    double c211 = a0 * c311 + a1 * c221 + a2 * c212;
    double d211d0 = c311;
    double d211d1 = c221;
    double d211d2 = c212;
     
    double c202 = a0 * c302 + a1 * c212 + a2 * c203;
    double d202d0 = c302;
    double d202d1 = c212;
    double d202d2 = c203;
    
    double c130 = a0 * c230 + a1 * c140 + a2 * c131;
    double d130d0 = c230;
    double d130d1 = c140;
    double d130d2 = c131;
    
    double c121 = a0 * c221 + a1 * c131 + a2 * c122;
    double d121d0 = c221;
    double d121d1 = c131;
    double d121d2 = c122;
    
    double c112 = a0 * c212 + a1 * c122 + a2 * c113;
    double d112d0 = c212;
    double d112d1 = c122;
    double d112d2 = c113;
     
    double c103 = a0 * c203 + a1 * c113 + a2 * c104;
    double d103d0 = c203;
    double d103d1 = c113;
    double d103d2 = c104;
    
    double c031 = a0 * c131 + a1 * c041 + a2 * c032;
    double d031d0 = c131;
    double d031d1 = c041;
    double d031d2 = c032;
    
    double c022 = a0 * c122 + a1 * c032 + a2 * c023;
    double d022d0 = c122;
    double d022d1 = c032;
    double d022d2 = c023;
    
    double c013 = a0 * c113 + a1 * c023 + a2 * c014;
    double d013d0 = c113;
    double d013d1 = c023;
    double d013d2 = c014;
    
    double c040 = a0 * c140 + a1 * c050 + a2 * c041;
    double d040d0 = c140;
    double d040d1 = c050;
    double d040d2 = c041;
    
    double c004 = a0 * c104 + a1 * c014 + a2 * c005;
    double d004d0 = c104;
    double d004d1 = c014;
    double d004d2 = c005;

    /* Layer 3 */
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double d300d0 = c400 + a0*d400d0 + a1*d310d0 + a2*d301d0;
    double d300d1 = c310 + a0*d400d1 + a1*d310d1 + a2*d301d1;
    double d300d2 = c301 + a0*d400d2 + a1*d310d2 + a2*d301d2;
    
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double d210d0 = c310 + a0*d310d0 + a1*d220d0 + a2*d211d0;
    double d210d1 = c220 + a0*d310d1 + a1*d220d1 + a2*d211d1;
    double d210d2 = c211 + a0*d310d2 + a1*d220d2 + a2*d211d2;
    
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double d201d0 = c301 + a0*d301d0 + a1*d211d0 + a2*d202d0;
    double d201d1 = c211 + a0*d301d1 + a1*d211d1 + a2*d202d1;
    double d201d2 = c202 + a0*d301d2 + a1*d211d2 + a2*d202d2;
    
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double d120d0 = c220 + a0*d220d0 + a1*d130d0 + a2*d121d0;
    double d120d1 = c130 + a0*d220d1 + a1*d130d1 + a2*d121d1;
    double d120d2 = c121 + a0*d220d2 + a1*d130d2 + a2*d121d2;
    
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double d111d0 = c211 + a0*d211d0 + a1*d121d0 + a2*d112d0;
    double d111d1 = c121 + a0*d211d1 + a1*d121d1 + a2*d112d1;
    double d111d2 = c112 + a0*d211d2 + a1*d121d2 + a2*d112d2;
    
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double d102d0 = c202 + a0*d202d0 + a1*d112d0 + a2*d103d0;
    double d102d1 = c112 + a0*d202d1 + a1*d112d1 + a2*d103d1;
    double d102d2 = c103 + a0*d202d2 + a1*d112d2 + a2*d103d2;
    
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double d030d0 = c130 + a0*d130d0 + a1*d040d0 + a2*d031d0;
    double d030d1 = c040 + a0*d130d1 + a1*d040d1 + a2*d031d1;
    double d030d2 = c031 + a0*d130d2 + a1*d040d2 + a2*d031d2;
    
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double d021d0 = c121 + a0*d121d0 + a1*d031d0 + a2*d022d0;
    double d021d1 = c031 + a0*d121d1 + a1*d031d1 + a2*d022d1;
    double d021d2 = c022 + a0*d121d2 + a1*d031d2 + a2*d022d2;
    
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double d012d0 = c112 + a0*d112d0 + a1*d022d0 + a2*d013d0;
    double d012d1 = c022 + a0*d112d1 + a1*d022d1 + a2*d013d1;
    double d012d2 = c013 + a0*d112d2 + a1*d022d2 + a2*d013d2;
    
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    double d003d0 = c103 + a0*d103d0 + a1*d013d0 + a2*d004d0;
    double d003d1 = c013 + a0*d103d1 + a1*d013d1 + a2*d004d1;
    double d003d2 = c004 + a0*d103d2 + a1*d013d2 + a2*d004d2;

    /* Layer 2 */
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double d200d0 = c300 + a0*d300d0 + a1*d210d0 + a2*d201d0;
    double d200d1 = c210 + a0*d300d1 + a1*d210d1 + a2*d201d1;
    double d200d2 = c201 + a0*d300d2 + a1*d210d2 + a2*d201d2;
    
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double d110d0 = c210 + a0*d210d0 + a1*d120d0 + a2*d111d0;
    double d110d1 = c120 + a0*d210d1 + a1*d120d1 + a2*d111d1;
    double d110d2 = c111 + a0*d210d2 + a1*d120d2 + a2*d111d2;
    
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double d101d0 = c201 + a0*d201d0 + a1*d111d0 + a2*d102d0;
    double d101d1 = c111 + a0*d201d1 + a1*d111d1 + a2*d102d1;
    double d101d2 = c102 + a0*d201d2 + a1*d111d2 + a2*d102d2;
    
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double d020d0 = c120 + a0*d120d0 + a1*d030d0 + a2*d021d0;
    double d020d1 = c030 + a0*d120d1 + a1*d030d1 + a2*d021d1;
    double d020d2 = c021 + a0*d120d2 + a1*d030d2 + a2*d021d2;
    
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double d011d0 = c111 + a0*d111d0 + a1*d021d0 + a2*d012d0;
    double d011d1 = c021 + a0*d111d1 + a1*d021d1 + a2*d012d1;
    double d011d2 = c012 + a0*d111d2 + a1*d021d2 + a2*d012d2;
    
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    double d002d0 = c102 + a0*d102d0 + a1*d012d0 + a2*d003d0;
    double d002d1 = c012 + a0*d102d1 + a1*d012d1 + a2*d003d1;
    double d002d2 = c003 + a0*d102d2 + a1*d012d2 + a2*d003d2;
    
    /* Layer 1 */
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double d100d0 = c200 + a0*d200d0 + a1*d110d0 + a2*d101d0;
    double d100d1 = c110 + a0*d200d1 + a1*d110d1 + a2*d101d1;
    double d100d2 = c101 + a0*d200d2 + a1*d110d2 + a2*d101d2;
          
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double d010d0 = c110 + a0*d110d0 + a1*d020d0 + a2*d011d0;
    double d010d1 = c020 + a0*d110d1 + a1*d020d1 + a2*d011d1;
    double d010d2 = c011 + a0*d110d2 + a1*d020d2 + a2*d011d2;
    
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    double d001d0 = c101 + a0*d101d0 + a1*d011d0 + a2*d002d0;
    double d001d1 = c011 + a0*d101d1 + a1*d011d1 + a2*d002d1;
    double d001d2 = c002 + a0*d101d2 + a1*d011d2 + a2*d002d2;

    /* Layer 0 */
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;
    double d000d0 = c100 + a0*d100d0 + a1*d010d0 + a2*d001d0;
    double d000d1 = c010 + a0*d100d1 + a1*d010d1 + a2*d001d1;
    double d000d2 = c001 + a0*d100d2 + a1*d010d2 + a2*d001d2;

    (*f) = c000;
    (*df) = (r3_t){{d000d0, d000d1, d000d2}};
  }
   
void SPDeCasteljau_Grad6(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df)
  {
    affirm(coeff.ne == 28, "wrong num of coeffs");
    
    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c600 = coeff.e[ 0];
    double c510 = coeff.e[ 1];
    double c501 = coeff.e[ 2];
    double c420 = coeff.e[ 3];
    double c411 = coeff.e[ 4];
    double c402 = coeff.e[ 5];
    double c330 = coeff.e[ 6];
    double c321 = coeff.e[ 7];
    double c312 = coeff.e[ 8];
    double c303 = coeff.e[ 9];
    double c240 = coeff.e[10];
    double c231 = coeff.e[11];
    double c222 = coeff.e[12];
    double c213 = coeff.e[13];
    double c204 = coeff.e[14];
    double c150 = coeff.e[15];
    double c141 = coeff.e[16];
    double c132 = coeff.e[17];
    double c123 = coeff.e[18];
    double c114 = coeff.e[19];
    double c105 = coeff.e[20];
    double c060 = coeff.e[21];
    double c051 = coeff.e[22];
    double c042 = coeff.e[23];
    double c033 = coeff.e[24];
    double c024 = coeff.e[25];
    double c015 = coeff.e[26];
    double c006 = coeff.e[27];
    
    
    /* dXYZdK is the derivative of cXYZ with respect to aK. */
    
    /* Layer 5 */
    
    double c500 = a0 * c600 + a1 * c510 + a2 * c501;
    double d500d0 = c600;
    double d500d1 = c510;
    double d500d2 = c501;
    
    double c410 = a0 * c510 + a1 * c420 + a2 * c411;
    double d410d0 = c510;
    double d410d1 = c420;
    double d410d2 = c411;
    
    double c401 = a0 * c501 + a1 * c411 + a2 * c402;
    double d401d0 = c501;
    double d401d1 = c411;
    double d401d2 = c402;
    
    double c320 = a0 * c420 + a1 * c330 + a2 * c321;
    double d320d0 = c420;
    double d320d1 = c330;
    double d320d2 = c321;
    
    double c311 = a0 * c411 + a1 * c321 + a2 * c312;
    double d311d0 = c411;
    double d311d1 = c321;
    double d311d2 = c312;
    
    double c302 = a0 * c402 + a1 * c312 + a2 * c303;
    double d302d0 = c402;
    double d302d1 = c312;
    double d302d2 = c303;
    
    double c230 = a0 * c330 + a1 * c240 + a2 * c231;
    double d230d0 = c330;
    double d230d1 = c240;
    double d230d2 = c231;
    
    double c221 = a0 * c321 + a1 * c231 + a2 * c222;
    double d221d0 = c321;
    double d221d1 = c231;
    double d221d2 = c222;
    
    double c212 = a0 * c312 + a1 * c222 + a2 * c213;
    double d212d0 = c312;
    double d212d1 = c222;
    double d212d2 = c213;
    
    double c203 = a0 * c303 + a1 * c213 + a2 * c204;
    double d203d0 = c303;
    double d203d1 = c213;
    double d203d2 = c204;
    
    double c140 = a0 * c240 + a1 * c150 + a2 * c141;
    double d140d0 = c240;
    double d140d1 = c150;
    double d140d2 = c141;
    
    double c131 = a0 * c231 + a1 * c141 + a2 * c132;
    double d131d0 = c231;
    double d131d1 = c141;
    double d131d2 = c132;
    
    double c122 = a0 * c222 + a1 * c132 + a2 * c123;
    double d122d0 = c222;
    double d122d1 = c132;
    double d122d2 = c123;
    
    double c113 = a0 * c213 + a1 * c123 + a2 * c114;
    double d113d0 = c213;
    double d113d1 = c123;
    double d113d2 = c114;
    
    double c104 = a0 * c204 + a1 * c114 + a2 * c105;
    double d104d0 = c204;
    double d104d1 = c114;
    double d104d2 = c105;
    
    double c050 = a0 * c150 + a1 * c060 + a2 * c051;
    double d050d0 = c150;
    double d050d1 = c060;
    double d050d2 = c051;
    
    double c041 = a0 * c141 + a1 * c051 + a2 * c042;
    double d041d0 = c141;
    double d041d1 = c051;
    double d041d2 = c042;
    
    double c032 = a0 * c132 + a1 * c042 + a2 * c033;
    double d032d0 = c132;
    double d032d1 = c042;
    double d032d2 = c033;
    
    double c023 = a0 * c123 + a1 * c033 + a2 * c024;
    double d023d0 = c123;
    double d023d1 = c033;
    double d023d2 = c024;
    
    double c014 = a0 * c114 + a1 * c024 + a2 * c015;
    double d014d0 = c114;
    double d014d1 = c024;
    double d014d2 = c015;
    
    double c005 = a0 * c105 + a1 * c015 + a2 * c006;
    double d005d0 = c105;
    double d005d1 = c015;
    double d005d2 = c006;
    
     /* Layer 4 */
    
    double c400 = a0 * c500 + a1 * c410 + a2 * c401;
    double d400d0 = c500 + a0*d500d0 + a1*d410d0 + a2*d401d0;
    double d400d1 = c410 + a0*d500d1 + a1*d410d1 + a2*d401d1;
    double d400d2 = c401 + a0*d500d2 + a1*d410d2 + a2*d401d2;
    
    double c310 = a0 * c410 + a1 * c320 + a2 * c311;
    double d310d0 = c410 + a0*d410d0 + a1*d320d0 + a2*d311d0;
    double d310d1 = c320 + a0*d410d1 + a1*d320d1 + a2*d311d1;
    double d310d2 = c311 + a0*d410d2 + a1*d320d2 + a2*d311d2;
    
    double c301 = a0 * c401 + a1 * c311 + a2 * c302;
    double d301d0 = c401 + a0*d401d0 + a1*d311d0 + a2*d302d0;
    double d301d1 = c311 + a0*d401d1 + a1*d311d1 + a2*d302d1;
    double d301d2 = c302 + a0*d401d2 + a1*d311d2 + a2*d302d2;
    
    double c220 = a0 * c320 + a1 * c230 + a2 * c221;
    double d220d0 = c320 + a0*d320d0 + a1*d230d0 + a2*d221d0;
    double d220d1 = c230 + a0*d320d1 + a1*d230d1 + a2*d221d1;
    double d220d2 = c221 + a0*d320d2 + a1*d230d2 + a2*d221d2;
    
    double c211 = a0 * c311 + a1 * c221 + a2 * c212;
    double d211d0 = c311 + a0*d311d0 + a1*d221d0 + a2*d212d0;
    double d211d1 = c221 + a0*d311d1 + a1*d221d1 + a2*d212d1;
    double d211d2 = c212 + a0*d311d2 + a1*d221d2 + a2*d212d2;
    
    double c202 = a0 * c302 + a1 * c212 + a2 * c203;
    double d202d0 = c302 + a0*d302d0 + a1*d212d0 + a2*d203d0;
    double d202d1 = c212 + a0*d302d1 + a1*d212d1 + a2*d203d1;
    double d202d2 = c203 + a0*d302d2 + a1*d212d2 + a2*d203d2;
    
    double c130 = a0 * c230 + a1 * c140 + a2 * c131;
    double d130d0 = c230 + a0*d230d0 + a1*d140d0 + a2*d131d0;
    double d130d1 = c140 + a0*d230d1 + a1*d140d1 + a2*d131d1;
    double d130d2 = c131 + a0*d230d2 + a1*d140d2 + a2*d131d2;
    
    double c121 = a0 * c221 + a1 * c131 + a2 * c122;
    double d121d0 = c221 + a0*d221d0 + a1*d131d0 + a2*d122d0;
    double d121d1 = c131 + a0*d221d1 + a1*d131d1 + a2*d122d1;
    double d121d2 = c122 + a0*d221d2 + a1*d131d2 + a2*d122d2;
    
    double c112 = a0 * c212 + a1 * c122 + a2 * c113;
    double d112d0 = c212 + a0*d212d0 + a1*d122d0 + a2*d113d0;
    double d112d1 = c122 + a0*d212d1 + a1*d122d1 + a2*d113d1;
    double d112d2 = c113 + a0*d212d2 + a1*d122d2 + a2*d113d2;
    
    double c103 = a0 * c203 + a1 * c113 + a2 * c104;
    double d103d0 = c203 + a0*d203d0 + a1*d113d0 + a2*d104d0;
    double d103d1 = c113 + a0*d203d1 + a1*d113d1 + a2*d104d1;
    double d103d2 = c104 + a0*d203d2 + a1*d113d2 + a2*d104d2;
    
    double c040 = a0 * c140 + a1 * c050 + a2 * c041;
    double d040d0 = c140 + a0*d140d0 + a1*d050d0 + a2*d041d0;
    double d040d1 = c050 + a0*d140d1 + a1*d050d1 + a2*d041d1;
    double d040d2 = c041 + a0*d140d2 + a1*d050d2 + a2*d041d2;
    
    double c031 = a0 * c131 + a1 * c041 + a2 * c032;
    double d031d0 = c131 + a0*d131d0 + a1*d041d0 + a2*d032d0;
    double d031d1 = c041 + a0*d131d1 + a1*d041d1 + a2*d032d1;
    double d031d2 = c032 + a0*d131d2 + a1*d041d2 + a2*d032d2;
    
    double c022 = a0 * c122 + a1 * c032 + a2 * c023;
    double d022d0 = c122 + a0*d122d0 + a1*d032d0 + a2*d023d0;
    double d022d1 = c032 + a0*d122d1 + a1*d032d1 + a2*d023d1;
    double d022d2 = c023 + a0*d122d2 + a1*d032d2 + a2*d023d2;
    
    double c013 = a0 * c113 + a1 * c023 + a2 * c014;
    double d013d0 = c113 + a0*d113d0 + a1*d023d0 + a2*d014d0;
    double d013d1 = c023 + a0*d113d1 + a1*d023d1 + a2*d014d1;
    double d013d2 = c014 + a0*d113d2 + a1*d023d2 + a2*d014d2;
    
    double c004 = a0 * c104 + a1 * c014 + a2 * c005;
    double d004d0 = c104 + a0*d104d0 + a1*d014d0 + a2*d005d0;
    double d004d1 = c014 + a0*d104d1 + a1*d014d1 + a2*d005d1;
    double d004d2 = c005 + a0*d104d2 + a1*d014d2 + a2*d005d2;
    
    
    /* Layer 3 */
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double d300d0 = c400 + a0*d400d0 + a1*d310d0 + a2*d301d0;
    double d300d1 = c310 + a0*d400d1 + a1*d310d1 + a2*d301d1;
    double d300d2 = c301 + a0*d400d2 + a1*d310d2 + a2*d301d2;
    
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double d210d0 = c310 + a0*d310d0 + a1*d220d0 + a2*d211d0;
    double d210d1 = c220 + a0*d310d1 + a1*d220d1 + a2*d211d1;
    double d210d2 = c211 + a0*d310d2 + a1*d220d2 + a2*d211d2;
    
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double d201d0 = c301 + a0*d301d0 + a1*d211d0 + a2*d202d0;
    double d201d1 = c211 + a0*d301d1 + a1*d211d1 + a2*d202d1;
    double d201d2 = c202 + a0*d301d2 + a1*d211d2 + a2*d202d2;
    
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double d120d0 = c220 + a0*d220d0 + a1*d130d0 + a2*d121d0;
    double d120d1 = c130 + a0*d220d1 + a1*d130d1 + a2*d121d1;
    double d120d2 = c121 + a0*d220d2 + a1*d130d2 + a2*d121d2;
    
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double d111d0 = c211 + a0*d211d0 + a1*d121d0 + a2*d112d0;
    double d111d1 = c121 + a0*d211d1 + a1*d121d1 + a2*d112d1;
    double d111d2 = c112 + a0*d211d2 + a1*d121d2 + a2*d112d2;
    
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double d102d0 = c202 + a0*d202d0 + a1*d112d0 + a2*d103d0;
    double d102d1 = c112 + a0*d202d1 + a1*d112d1 + a2*d103d1;
    double d102d2 = c103 + a0*d202d2 + a1*d112d2 + a2*d103d2;
    
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double d030d0 = c130 + a0*d130d0 + a1*d040d0 + a2*d031d0;
    double d030d1 = c040 + a0*d130d1 + a1*d040d1 + a2*d031d1;
    double d030d2 = c031 + a0*d130d2 + a1*d040d2 + a2*d031d2;
    
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double d021d0 = c121 + a0*d121d0 + a1*d031d0 + a2*d022d0;
    double d021d1 = c031 + a0*d121d1 + a1*d031d1 + a2*d022d1;
    double d021d2 = c022 + a0*d121d2 + a1*d031d2 + a2*d022d2;
    
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double d012d0 = c112 + a0*d112d0 + a1*d022d0 + a2*d013d0;
    double d012d1 = c022 + a0*d112d1 + a1*d022d1 + a2*d013d1;
    double d012d2 = c013 + a0*d112d2 + a1*d022d2 + a2*d013d2;
    
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    double d003d0 = c103 + a0*d103d0 + a1*d013d0 + a2*d004d0;
    double d003d1 = c013 + a0*d103d1 + a1*d013d1 + a2*d004d1;
    double d003d2 = c004 + a0*d103d2 + a1*d013d2 + a2*d004d2;
     
    /* Layer 2 */
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double d200d0 = c300 + a0*d300d0 + a1*d210d0 + a2*d201d0;
    double d200d1 = c210 + a0*d300d1 + a1*d210d1 + a2*d201d1;
    double d200d2 = c201 + a0*d300d2 + a1*d210d2 + a2*d201d2;
    
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double d110d0 = c210 + a0*d210d0 + a1*d120d0 + a2*d111d0;
    double d110d1 = c120 + a0*d210d1 + a1*d120d1 + a2*d111d1;
    double d110d2 = c111 + a0*d210d2 + a1*d120d2 + a2*d111d2;
    
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double d101d0 = c201 + a0*d201d0 + a1*d111d0 + a2*d102d0;
    double d101d1 = c111 + a0*d201d1 + a1*d111d1 + a2*d102d1;
    double d101d2 = c102 + a0*d201d2 + a1*d111d2 + a2*d102d2;
    
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double d020d0 = c120 + a0*d120d0 + a1*d030d0 + a2*d021d0;
    double d020d1 = c030 + a0*d120d1 + a1*d030d1 + a2*d021d1;
    double d020d2 = c021 + a0*d120d2 + a1*d030d2 + a2*d021d2;
    
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double d011d0 = c111 + a0*d111d0 + a1*d021d0 + a2*d012d0;
    double d011d1 = c021 + a0*d111d1 + a1*d021d1 + a2*d012d1;
    double d011d2 = c012 + a0*d111d2 + a1*d021d2 + a2*d012d2;
    
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    double d002d0 = c102 + a0*d102d0 + a1*d012d0 + a2*d003d0;
    double d002d1 = c012 + a0*d102d1 + a1*d012d1 + a2*d003d1;
    double d002d2 = c003 + a0*d102d2 + a1*d012d2 + a2*d003d2;
    
    /* Layer 1 */
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double d100d0 = c200 + a0*d200d0 + a1*d110d0 + a2*d101d0;
    double d100d1 = c110 + a0*d200d1 + a1*d110d1 + a2*d101d1;
    double d100d2 = c101 + a0*d200d2 + a1*d110d2 + a2*d101d2;
          
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double d010d0 = c110 + a0*d110d0 + a1*d020d0 + a2*d011d0;
    double d010d1 = c020 + a0*d110d1 + a1*d020d1 + a2*d011d1;
    double d010d2 = c011 + a0*d110d2 + a1*d020d2 + a2*d011d2;
    
    
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    double d001d0 = c101 + a0*d101d0 + a1*d011d0 + a2*d002d0;
    double d001d1 = c011 + a0*d101d1 + a1*d011d1 + a2*d002d1;
    double d001d2 = c002 + a0*d101d2 + a1*d011d2 + a2*d002d2;
    
    /* Layer 0 */
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;
    double d000d0 = c100 + a0*d100d0 + a1*d010d0 + a2*d001d0;
    double d000d1 = c010 + a0*d100d1 + a1*d010d1 + a2*d001d1;
    double d000d2 = c001 + a0*d100d2 + a1*d010d2 + a2*d001d2;

    (*f) = c000;
    (*df) = (r3_t){{d000d0, d000d1, d000d2}};
  }

void SPDeCasteljau_Grad7(BezCoeff_vec_t coeff, r3_t *a, double *f, r3_t *df)
  {
    affirm(coeff.ne == 36, "wrong num of coeffs");

    double a0 = a->c[0];
    double a1 = a->c[1];
    double a2 = a->c[2];
    
    double c700 = coeff.e[ 0];
    double c610 = coeff.e[ 1];
    double c601 = coeff.e[ 2];
    double c520 = coeff.e[ 3];
    double c511 = coeff.e[ 4];
    double c502 = coeff.e[ 5];
    double c430 = coeff.e[ 6];
    double c421 = coeff.e[ 7];
    double c412 = coeff.e[ 8];
    double c403 = coeff.e[ 9];
    double c340 = coeff.e[10];
    double c331 = coeff.e[11];
    double c322 = coeff.e[12];
    double c313 = coeff.e[13];
    double c304 = coeff.e[14];
    double c250 = coeff.e[15];
    double c241 = coeff.e[16];
    double c232 = coeff.e[17];
    double c223 = coeff.e[18];
    double c214 = coeff.e[19];
    double c205 = coeff.e[20];
    double c160 = coeff.e[21];
    double c151 = coeff.e[22];
    double c142 = coeff.e[23];
    double c133 = coeff.e[24];
    double c124 = coeff.e[25];
    double c115 = coeff.e[26];
    double c106 = coeff.e[27];
    double c070 = coeff.e[28];
    double c061 = coeff.e[29];
    double c052 = coeff.e[30];
    double c043 = coeff.e[31];
    double c034 = coeff.e[32];
    double c025 = coeff.e[33];
    double c016 = coeff.e[34];
    double c007 = coeff.e[35];

    /* dXYZdK is the derivative of cXYZ with respect to aK. */
    
    /* Layer 6 */
    
    double c600 = a0 * c700 + a1 * c610 + a2 * c601;
    double d600d0 = c700;
    double d600d1 = c610;
    double d600d2 = c601;

    double c510 = a0 * c610 + a1 * c520 + a2 * c511;
    double d510d0 = c610;
    double d510d1 = c520;
    double d510d2 = c511;
    
    double c501 = a0 * c601 + a1 * c511 + a2 * c502;
    double d501d0 = c601;
    double d501d1 = c511;
    double d501d2 = c502;
    
    double c420 = a0 * c520 + a1 * c430 + a2 * c421;
    double d420d0 = c520;
    double d420d1 = c430;
    double d420d2 = c421;
    
    double c411 = a0 * c511 + a1 * c421 + a2 * c412;
    double d411d0 = c511;
    double d411d1 = c421;
    double d411d2 = c412;
    
    double c402 = a0 * c502 + a1 * c412 + a2 * c403;
    double d402d0 = c502;
    double d402d1 = c412;
    double d402d2 = c403;
    
    double c330 = a0 * c430 + a1 * c340 + a2 * c331;
    double d330d0 = c430;
    double d330d1 = c340;
    double d330d2 = c331;
    
    double c321 = a0 * c421 + a1 * c331 + a2 * c322;
    double d321d0 = c421;
    double d321d1 = c331;
    double d321d2 = c322;
    
    double c312 = a0 * c412 + a1 * c322 + a2 * c313;
    double d312d0 = c412;
    double d312d1 = c322;
    double d312d2 = c313;
    
    double c303 = a0 * c403 + a1 * c313 + a2 * c304;
    double d303d0 = c403;
    double d303d1 = c313;
    double d303d2 = c304;
    
    double c240 = a0 * c340 + a1 * c250 + a2 * c241;
    double d240d0 = c340;
    double d240d1 = c250;
    double d240d2 = c241;
    
    double c231 = a0 * c331 + a1 * c241 + a2 * c232;
    double d231d0 = c331;
    double d231d1 = c241;
    double d231d2 = c232;
    
    double c222 = a0 * c322 + a1 * c232 + a2 * c223;
    double d222d0 = c322;
    double d222d1 = c232;
    double d222d2 = c223;
    
    double c213 = a0 * c313 + a1 * c223 + a2 * c214;
    double d213d0 = c313;
    double d213d1 = c223;
    double d213d2 = c214;
    
    double c204 = a0 * c304 + a1 * c214 + a2 * c205;
    double d204d0 = c304;
    double d204d1 = c214;
    double d204d2 = c205;
    
    double c150 = a0 * c250 + a1 * c160 + a2 * c151;
    double d150d0 = c250;
    double d150d1 = c160;
    double d150d2 = c151;
    
    double c141 = a0 * c241 + a1 * c151 + a2 * c142;
    double d141d0 = c241;
    double d141d1 = c151;
    double d141d2 = c142;
    
    double c132 = a0 * c232 + a1 * c142 + a2 * c133;
    double d132d0 = c232;
    double d132d1 = c142;
    double d132d2 = c133;
    
    double c123 = a0 * c223 + a1 * c133 + a2 * c124;
    double d123d0 = c223;
    double d123d1 = c133;
    double d123d2 = c124;
    
    double c114 = a0 * c214 + a1 * c124 + a2 * c115;
    double d114d0 = c214;
    double d114d1 = c124;
    double d114d2 = c115;
    
    double c105 = a0 * c205 + a1 * c115 + a2 * c106;
    double d105d0 = c205;
    double d105d1 = c115;
    double d105d2 = c106;
    
    double c060 = a0 * c160 + a1 * c070 + a2 * c061;
    double d060d0 = c160;
    double d060d1 = c070;
    double d060d2 = c061;
    
    double c051 = a0 * c151 + a1 * c061 + a2 * c052;
    double d051d0 = c151;
    double d051d1 = c061;
    double d051d2 = c052;
    
    double c042 = a0 * c142 + a1 * c052 + a2 * c043;
    double d042d0 = c142;
    double d042d1 = c052;
    double d042d2 = c043;
    
    double c033 = a0 * c133 + a1 * c043 + a2 * c034;
    double d033d0 = c133;
    double d033d1 = c043;
    double d033d2 = c034;
    
    double c024 = a0 * c124 + a1 * c034 + a2 * c025;
    double d024d0 = c124;
    double d024d1 = c034;
    double d024d2 = c025;
    
    double c015 = a0 * c115 + a1 * c025 + a2 * c016;
    double d015d0 = c115;
    double d015d1 = c025;
    double d015d2 = c016;
    
    double c006 = a0 * c106 + a1 * c016 + a2 * c007;
    double d006d0 = c106;
    double d006d1 = c016;
    double d006d2 = c007;
    
    /* Layer 5 */
    
    double c500 = a0 * c600 + a1 * c510 + a2 * c501;
    double d500d0 = c600 + a0*d600d0 + a1*d510d0 + a2*d501d0;
    double d500d1 = c510 + a0*d600d1 + a1*d510d1 + a2*d501d1;
    double d500d2 = c501 + a0*d600d2 + a1*d510d2 + a2*d501d2;
    
    double c410 = a0 * c510 + a1 * c420 + a2 * c411;
    double d410d0 = c510 + a0*d510d0 + a1*d420d0 + a2*d411d0;
    double d410d1 = c420 + a0*d510d1 + a1*d420d1 + a2*d411d1;
    double d410d2 = c411 + a0*d510d2 + a1*d420d2 + a2*d411d2;
    
    double c401 = a0 * c501 + a1 * c411 + a2 * c402;
    double d401d0 = c501 + a0*d501d0 + a1*d411d0 + a2*d402d0;
    double d401d1 = c411 + a0*d501d1 + a1*d411d1 + a2*d402d1;
    double d401d2 = c402 + a0*d501d2 + a1*d411d2 + a2*d402d2;
    
    double c320 = a0 * c420 + a1 * c330 + a2 * c321;
    double d320d0 = c420 + a0*d420d0 + a1*d330d0 + a2*d321d0;
    double d320d1 = c330 + a0*d420d1 + a1*d330d1 + a2*d321d1;
    double d320d2 = c321 + a0*d420d2 + a1*d330d2 + a2*d321d2;
    
    double c311 = a0 * c411 + a1 * c321 + a2 * c312;
    double d311d0 = c411 + a0*d411d0 + a1*d321d0 + a2*d312d0;
    double d311d1 = c321 + a0*d411d1 + a1*d321d1 + a2*d312d1;
    double d311d2 = c312 + a0*d411d2 + a1*d321d2 + a2*d312d2;
    
    double c302 = a0 * c402 + a1 * c312 + a2 * c303;
    double d302d0 = c402 + a0*d402d0 + a1*d312d0 + a2*d303d0;
    double d302d1 = c312 + a0*d402d1 + a1*d312d1 + a2*d303d1;
    double d302d2 = c303 + a0*d402d2 + a1*d312d2 + a2*d303d2;
    
    double c230 = a0 * c330 + a1 * c240 + a2 * c231;
    double d230d0 = c330 + a0*d330d0 + a1*d240d0 + a2*d231d0;
    double d230d1 = c240 + a0*d330d1 + a1*d240d1 + a2*d231d1;
    double d230d2 = c231 + a0*d330d2 + a1*d240d2 + a2*d231d2;
    
    double c221 = a0 * c321 + a1 * c231 + a2 * c222;
    double d221d0 = c321 + a0*d321d0 + a1*d231d0 + a2*d222d0;
    double d221d1 = c231 + a0*d321d1 + a1*d231d1 + a2*d222d1;
    double d221d2 = c222 + a0*d321d2 + a1*d231d2 + a2*d222d2;
    
    double c212 = a0 * c312 + a1 * c222 + a2 * c213;
    double d212d0 = c312 + a0*d312d0 + a1*d222d0 + a2*d213d0;
    double d212d1 = c222 + a0*d312d1 + a1*d222d1 + a2*d213d1;
    double d212d2 = c213 + a0*d312d2 + a1*d222d2 + a2*d213d2;
    
    double c203 = a0 * c303 + a1 * c213 + a2 * c204;
    double d203d0 = c303 + a0*d303d0 + a1*d213d0 + a2*d204d0;
    double d203d1 = c213 + a0*d303d1 + a1*d213d1 + a2*d204d1;
    double d203d2 = c204 + a0*d303d2 + a1*d213d2 + a2*d204d2;
    
    double c140 = a0 * c240 + a1 * c150 + a2 * c141;
    double d140d0 = c240 + a0*d240d0 + a1*d150d0 + a2*d141d0;
    double d140d1 = c150 + a0*d240d1 + a1*d150d1 + a2*d141d1;
    double d140d2 = c141 + a0*d240d2 + a1*d150d2 + a2*d141d2;
    
    double c131 = a0 * c231 + a1 * c141 + a2 * c132;
    double d131d0 = c231 + a0*d231d0 + a1*d141d0 + a2*d132d0;
    double d131d1 = c141 + a0*d231d1 + a1*d141d1 + a2*d132d1;
    double d131d2 = c132 + a0*d231d2 + a1*d141d2 + a2*d132d2;
    
    double c122 = a0 * c222 + a1 * c132 + a2 * c123;
    double d122d0 = c222 + a0*d222d0 + a1*d132d0 + a2*d123d0;
    double d122d1 = c132 + a0*d222d1 + a1*d132d1 + a2*d123d1;
    double d122d2 = c123 + a0*d222d2 + a1*d132d2 + a2*d123d2;
    
    double c113 = a0 * c213 + a1 * c123 + a2 * c114;
    double d113d0 = c213 + a0*d213d0 + a1*d123d0 + a2*d114d0;
    double d113d1 = c123 + a0*d213d1 + a1*d123d1 + a2*d114d1;
    double d113d2 = c114 + a0*d213d2 + a1*d123d2 + a2*d114d2;
    
    double c104 = a0 * c204 + a1 * c114 + a2 * c105;
    double d104d0 = c204 + a0*d204d0 + a1*d114d0 + a2*d105d0;
    double d104d1 = c114 + a0*d204d1 + a1*d114d1 + a2*d105d1;
    double d104d2 = c105 + a0*d204d2 + a1*d114d2 + a2*d105d2;
    
    double c050 = a0 * c150 + a1 * c060 + a2 * c051;
    double d050d0 = c150 + a0*d150d0 + a1*d060d0 + a2*d051d0;
    double d050d1 = c060 + a0*d150d1 + a1*d060d1 + a2*d051d1;
    double d050d2 = c051 + a0*d150d2 + a1*d060d2 + a2*d051d2;
    
    double c041 = a0 * c141 + a1 * c051 + a2 * c042;
    double d041d0 = c141 + a0*d141d0 + a1*d051d0 + a2*d042d0;
    double d041d1 = c051 + a0*d141d1 + a1*d051d1 + a2*d042d1;
    double d041d2 = c042 + a0*d141d2 + a1*d051d2 + a2*d042d2;
    
    double c032 = a0 * c132 + a1 * c042 + a2 * c033;
    double d032d0 = c132 + a0*d132d0 + a1*d042d0 + a2*d033d0;
    double d032d1 = c042 + a0*d132d1 + a1*d042d1 + a2*d033d1;
    double d032d2 = c033 + a0*d132d2 + a1*d042d2 + a2*d033d2;
    
    double c023 = a0 * c123 + a1 * c033 + a2 * c024;
    double d023d0 = c123 + a0*d123d0 + a1*d033d0 + a2*d024d0;
    double d023d1 = c033 + a0*d123d1 + a1*d033d1 + a2*d024d1;
    double d023d2 = c024 + a0*d123d2 + a1*d033d2 + a2*d024d2;
    
    double c014 = a0 * c114 + a1 * c024 + a2 * c015;
    double d014d0 = c114 + a0*d114d0 + a1*d024d0 + a2*d015d0;
    double d014d1 = c024 + a0*d114d1 + a1*d024d1 + a2*d015d1;
    double d014d2 = c015 + a0*d114d2 + a1*d024d2 + a2*d015d2;
    
    double c005 = a0 * c105 + a1 * c015 + a2 * c006;
    double d005d0 = c105 + a0*d105d0 + a1*d015d0 + a2*d006d0;
    double d005d1 = c015 + a0*d105d1 + a1*d015d1 + a2*d006d1;
    double d005d2 = c006 + a0*d105d2 + a1*d015d2 + a2*d006d2;
    
     /* Layer 4 */
    
    double c400 = a0 * c500 + a1 * c410 + a2 * c401;
    double d400d0 = c500 + a0*d500d0 + a1*d410d0 + a2*d401d0;
    double d400d1 = c410 + a0*d500d1 + a1*d410d1 + a2*d401d1;
    double d400d2 = c401 + a0*d500d2 + a1*d410d2 + a2*d401d2;
    
    double c310 = a0 * c410 + a1 * c320 + a2 * c311;
    double d310d0 = c410 + a0*d410d0 + a1*d320d0 + a2*d311d0;
    double d310d1 = c320 + a0*d410d1 + a1*d320d1 + a2*d311d1;
    double d310d2 = c311 + a0*d410d2 + a1*d320d2 + a2*d311d2;
    
    double c301 = a0 * c401 + a1 * c311 + a2 * c302;
    double d301d0 = c401 + a0*d401d0 + a1*d311d0 + a2*d302d0;
    double d301d1 = c311 + a0*d401d1 + a1*d311d1 + a2*d302d1;
    double d301d2 = c302 + a0*d401d2 + a1*d311d2 + a2*d302d2;
    
    double c220 = a0 * c320 + a1 * c230 + a2 * c221;
    double d220d0 = c320 + a0*d320d0 + a1*d230d0 + a2*d221d0;
    double d220d1 = c230 + a0*d320d1 + a1*d230d1 + a2*d221d1;
    double d220d2 = c221 + a0*d320d2 + a1*d230d2 + a2*d221d2;
    
    double c211 = a0 * c311 + a1 * c221 + a2 * c212;
    double d211d0 = c311 + a0*d311d0 + a1*d221d0 + a2*d212d0;
    double d211d1 = c221 + a0*d311d1 + a1*d221d1 + a2*d212d1;
    double d211d2 = c212 + a0*d311d2 + a1*d221d2 + a2*d212d2;
    
    double c202 = a0 * c302 + a1 * c212 + a2 * c203;
    double d202d0 = c302 + a0*d302d0 + a1*d212d0 + a2*d203d0;
    double d202d1 = c212 + a0*d302d1 + a1*d212d1 + a2*d203d1;
    double d202d2 = c203 + a0*d302d2 + a1*d212d2 + a2*d203d2;
    
    double c130 = a0 * c230 + a1 * c140 + a2 * c131;
    double d130d0 = c230 + a0*d230d0 + a1*d140d0 + a2*d131d0;
    double d130d1 = c140 + a0*d230d1 + a1*d140d1 + a2*d131d1;
    double d130d2 = c131 + a0*d230d2 + a1*d140d2 + a2*d131d2;
    
    double c121 = a0 * c221 + a1 * c131 + a2 * c122;
    double d121d0 = c221 + a0*d221d0 + a1*d131d0 + a2*d122d0;
    double d121d1 = c131 + a0*d221d1 + a1*d131d1 + a2*d122d1;
    double d121d2 = c122 + a0*d221d2 + a1*d131d2 + a2*d122d2;
    
    double c112 = a0 * c212 + a1 * c122 + a2 * c113;
    double d112d0 = c212 + a0*d212d0 + a1*d122d0 + a2*d113d0;
    double d112d1 = c122 + a0*d212d1 + a1*d122d1 + a2*d113d1;
    double d112d2 = c113 + a0*d212d2 + a1*d122d2 + a2*d113d2;
    
    double c103 = a0 * c203 + a1 * c113 + a2 * c104;
    double d103d0 = c203 + a0*d203d0 + a1*d113d0 + a2*d104d0;
    double d103d1 = c113 + a0*d203d1 + a1*d113d1 + a2*d104d1;
    double d103d2 = c104 + a0*d203d2 + a1*d113d2 + a2*d104d2;
    
    double c040 = a0 * c140 + a1 * c050 + a2 * c041;
    double d040d0 = c140 + a0*d140d0 + a1*d050d0 + a2*d041d0;
    double d040d1 = c050 + a0*d140d1 + a1*d050d1 + a2*d041d1;
    double d040d2 = c041 + a0*d140d2 + a1*d050d2 + a2*d041d2;
    
    double c031 = a0 * c131 + a1 * c041 + a2 * c032;
    double d031d0 = c131 + a0*d131d0 + a1*d041d0 + a2*d032d0;
    double d031d1 = c041 + a0*d131d1 + a1*d041d1 + a2*d032d1;
    double d031d2 = c032 + a0*d131d2 + a1*d041d2 + a2*d032d2;
    
    double c022 = a0 * c122 + a1 * c032 + a2 * c023;
    double d022d0 = c122 + a0*d122d0 + a1*d032d0 + a2*d023d0;
    double d022d1 = c032 + a0*d122d1 + a1*d032d1 + a2*d023d1;
    double d022d2 = c023 + a0*d122d2 + a1*d032d2 + a2*d023d2;
    
    double c013 = a0 * c113 + a1 * c023 + a2 * c014;
    double d013d0 = c113 + a0*d113d0 + a1*d023d0 + a2*d014d0;
    double d013d1 = c023 + a0*d113d1 + a1*d023d1 + a2*d014d1;
    double d013d2 = c014 + a0*d113d2 + a1*d023d2 + a2*d014d2;
    
    double c004 = a0 * c104 + a1 * c014 + a2 * c005;
    double d004d0 = c104 + a0*d104d0 + a1*d014d0 + a2*d005d0;
    double d004d1 = c014 + a0*d104d1 + a1*d014d1 + a2*d005d1;
    double d004d2 = c005 + a0*d104d2 + a1*d014d2 + a2*d005d2;
    
    
    /* Layer 3 */
    
    double c300 = a0 * c400 + a1 * c310 + a2 * c301;
    double d300d0 = c400 + a0*d400d0 + a1*d310d0 + a2*d301d0;
    double d300d1 = c310 + a0*d400d1 + a1*d310d1 + a2*d301d1;
    double d300d2 = c301 + a0*d400d2 + a1*d310d2 + a2*d301d2;
    
    double c210 = a0 * c310 + a1 * c220 + a2 * c211;
    double d210d0 = c310 + a0*d310d0 + a1*d220d0 + a2*d211d0;
    double d210d1 = c220 + a0*d310d1 + a1*d220d1 + a2*d211d1;
    double d210d2 = c211 + a0*d310d2 + a1*d220d2 + a2*d211d2;
    
    double c201 = a0 * c301 + a1 * c211 + a2 * c202;
    double d201d0 = c301 + a0*d301d0 + a1*d211d0 + a2*d202d0;
    double d201d1 = c211 + a0*d301d1 + a1*d211d1 + a2*d202d1;
    double d201d2 = c202 + a0*d301d2 + a1*d211d2 + a2*d202d2;
    
    double c120 = a0 * c220 + a1 * c130 + a2 * c121;
    double d120d0 = c220 + a0*d220d0 + a1*d130d0 + a2*d121d0;
    double d120d1 = c130 + a0*d220d1 + a1*d130d1 + a2*d121d1;
    double d120d2 = c121 + a0*d220d2 + a1*d130d2 + a2*d121d2;
    
    double c111 = a0 * c211 + a1 * c121 + a2 * c112;
    double d111d0 = c211 + a0*d211d0 + a1*d121d0 + a2*d112d0;
    double d111d1 = c121 + a0*d211d1 + a1*d121d1 + a2*d112d1;
    double d111d2 = c112 + a0*d211d2 + a1*d121d2 + a2*d112d2;
    
    double c102 = a0 * c202 + a1 * c112 + a2 * c103;
    double d102d0 = c202 + a0*d202d0 + a1*d112d0 + a2*d103d0;
    double d102d1 = c112 + a0*d202d1 + a1*d112d1 + a2*d103d1;
    double d102d2 = c103 + a0*d202d2 + a1*d112d2 + a2*d103d2;
    
    double c030 = a0 * c130 + a1 * c040 + a2 * c031;
    double d030d0 = c130 + a0*d130d0 + a1*d040d0 + a2*d031d0;
    double d030d1 = c040 + a0*d130d1 + a1*d040d1 + a2*d031d1;
    double d030d2 = c031 + a0*d130d2 + a1*d040d2 + a2*d031d2;
    
    double c021 = a0 * c121 + a1 * c031 + a2 * c022;
    double d021d0 = c121 + a0*d121d0 + a1*d031d0 + a2*d022d0;
    double d021d1 = c031 + a0*d121d1 + a1*d031d1 + a2*d022d1;
    double d021d2 = c022 + a0*d121d2 + a1*d031d2 + a2*d022d2;
    
    double c012 = a0 * c112 + a1 * c022 + a2 * c013;
    double d012d0 = c112 + a0*d112d0 + a1*d022d0 + a2*d013d0;
    double d012d1 = c022 + a0*d112d1 + a1*d022d1 + a2*d013d1;
    double d012d2 = c013 + a0*d112d2 + a1*d022d2 + a2*d013d2;
    
    double c003 = a0 * c103 + a1 * c013 + a2 * c004;
    double d003d0 = c103 + a0*d103d0 + a1*d013d0 + a2*d004d0;
    double d003d1 = c013 + a0*d103d1 + a1*d013d1 + a2*d004d1;
    double d003d2 = c004 + a0*d103d2 + a1*d013d2 + a2*d004d2;
     
    /* Layer 2 */
    
    double c200 = a0 * c300 + a1 * c210 + a2 * c201;
    double d200d0 = c300 + a0*d300d0 + a1*d210d0 + a2*d201d0;
    double d200d1 = c210 + a0*d300d1 + a1*d210d1 + a2*d201d1;
    double d200d2 = c201 + a0*d300d2 + a1*d210d2 + a2*d201d2;
    
    double c110 = a0 * c210 + a1 * c120 + a2 * c111;
    double d110d0 = c210 + a0*d210d0 + a1*d120d0 + a2*d111d0;
    double d110d1 = c120 + a0*d210d1 + a1*d120d1 + a2*d111d1;
    double d110d2 = c111 + a0*d210d2 + a1*d120d2 + a2*d111d2;
    
    double c101 = a0 * c201 + a1 * c111 + a2 * c102;
    double d101d0 = c201 + a0*d201d0 + a1*d111d0 + a2*d102d0;
    double d101d1 = c111 + a0*d201d1 + a1*d111d1 + a2*d102d1;
    double d101d2 = c102 + a0*d201d2 + a1*d111d2 + a2*d102d2;
    
    double c020 = a0 * c120 + a1 * c030 + a2 * c021;
    double d020d0 = c120 + a0*d120d0 + a1*d030d0 + a2*d021d0;
    double d020d1 = c030 + a0*d120d1 + a1*d030d1 + a2*d021d1;
    double d020d2 = c021 + a0*d120d2 + a1*d030d2 + a2*d021d2;
    
    double c011 = a0 * c111 + a1 * c021 + a2 * c012;
    double d011d0 = c111 + a0*d111d0 + a1*d021d0 + a2*d012d0;
    double d011d1 = c021 + a0*d111d1 + a1*d021d1 + a2*d012d1;
    double d011d2 = c012 + a0*d111d2 + a1*d021d2 + a2*d012d2;
    
    double c002 = a0 * c102 + a1 * c012 + a2 * c003;
    double d002d0 = c102 + a0*d102d0 + a1*d012d0 + a2*d003d0;
    double d002d1 = c012 + a0*d102d1 + a1*d012d1 + a2*d003d1;
    double d002d2 = c003 + a0*d102d2 + a1*d012d2 + a2*d003d2;
    
    /* Layer 1 */
    
    double c100 = a0 * c200 + a1 * c110 + a2 * c101;
    double d100d0 = c200 + a0*d200d0 + a1*d110d0 + a2*d101d0;
    double d100d1 = c110 + a0*d200d1 + a1*d110d1 + a2*d101d1;
    double d100d2 = c101 + a0*d200d2 + a1*d110d2 + a2*d101d2;
          
    double c010 = a0 * c110 + a1 * c020 + a2 * c011;
    double d010d0 = c110 + a0*d110d0 + a1*d020d0 + a2*d011d0;
    double d010d1 = c020 + a0*d110d1 + a1*d020d1 + a2*d011d1;
    double d010d2 = c011 + a0*d110d2 + a1*d020d2 + a2*d011d2;
    
    double c001 = a0 * c101 + a1 * c011 + a2 * c002;
    double d001d0 = c101 + a0*d101d0 + a1*d011d0 + a2*d002d0;
    double d001d1 = c011 + a0*d101d1 + a1*d011d1 + a2*d002d1;
    double d001d2 = c002 + a0*d101d2 + a1*d011d2 + a2*d002d2;
    
    /* Layer 0 */
    
    double c000 = a0 * c100 + a1 * c010 + a2 * c001;
    double d000d0 = c100 + a0*d100d0 + a1*d010d0 + a2*d001d0;
    double d000d1 = c010 + a0*d100d1 + a1*d010d1 + a2*d001d1;
    double d000d2 = c001 + a0*d100d2 + a1*d010d2 + a2*d001d2;

    (*f) = c000;
    (*df) = (r3_t){{d000d0, d000d1, d000d2}};
  }

