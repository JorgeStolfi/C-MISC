#define MAX_CARD 262144

typedef int Vertex;
typedef unsigned short int VertexSmall;

typedef struct enode {
  struct enode *next;
  struct enode *mate;
  Vertex other;
  Boolean matched;
} Edge;

typedef struct {
  Edge *current;
  Edge *e;
  Vertex matched; /* used as Boolean until check; then converted to Vertex */
  int excess;
  int label;
} Vvertex;

typedef struct {
  Edge *current;
  Edge *e;
  Vertex matchedTo;
  int label;
} Uvertex;

typedef struct {
  Uvertex A[MAX_CARD];
  Vertex card;
  int nEdges;
} Upartition;

typedef struct {
  Vvertex A[MAX_CARD];
  Vertex card;
} Vpartition;
