
#ifndef SystemIO_H
#define SystemIO_H



#include <Rd.h>
#include <Wr.h>
#include <r3.h>
#include <r3x3.h>

typedef
  Tetrahedron == struct ?? {
      nat p0, p1, p2, p3;
      r3x3_t A;
      double density, alpha, beta, eta1, eta2;
    }
  Face == struct ?? { u, v, w, nat tx, mu1, mu2, e: double; }
  Constraint == struct ?? { a: nat; }
  Vectors3D == r3_t_vec;
  PlainVectors3D == double_vec;
  
void ReadVectors(FILE *rd, nat n, pos, Vectors3D vel);
void ReadPlainVectors(FILE *rd, nat n, pos, PlainVectors3D vel);
void ReadFace(FILE *rd, VAR Face f);
void ReadFixed(FILE *rd, VAR nat i);
PROCEDURE ReadKinetic(FILE *rd, 
  double *?ta, tb; VAR nat k, VAR p, v, double xi,
);
void ReadTexture(FILE *rd, VAR Texture p);
PROCEDURE WriteVectors(FILE *wr, nat n, pos, Vectors3D vel, 
                       int base);
PROCEDURE WritePlainVectors(FILE *wr, nat n, pos, PlainVectors3D vel,
                            int base);
void WriteFace(FILE *wr, READONLY Face f, nat base);
void WriteFixed(FILE *wr, nat i, nat base);
void WriteKinetic(FILE *wr, ta, double tb, nat k, p, v, double xi);
void WriteTexture(FILE *wr, READONLY Texture p);
;
} /* SystemIO */.

#endif
