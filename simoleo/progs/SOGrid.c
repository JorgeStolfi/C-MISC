/* See SOGrid.h */
/* Last edited on 2005-06-05 20:58:45 by stolfi */

#include <SOGrid.h>

#include <dg_grid.h>

#include <affirm.h>
#include <nat.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <stdlib.h>

/* INTERNAL PROTOTYPES */

void SOGrid_subtree_write(FILE *wr, dg_tree_node_t *p);
  /* Writes the subtree rooted at {p} to {wr}. If {p} is 
    NULL, writes "*"; else writes "(A,B)" where "A" and "B"
    are the low and high subtrees of {p}. */

dg_tree_node_t *SOGrid_subtree_read(FILE *rd, dg_tree_node_t *p, dg_cell_index_t k);
  /* Reads a subtree from {rd}. If the subtree is not NULL, its root
    will have parent {p} and index {k}. */

/* IMPLEMENTATIONS */

/* ALLOCATION */

SOGrid_Tree *SOGrid_Tree_new(dg_dim_t d)
  { void *v = notnull(malloc(sizeof(SOGrid_Tree)), "no mem for SOGrid_Tree");
    SOGrid_Tree *h = (SOGrid_Tree*)v;
    h->d = d;
    h->root = dg_tree_node_new(NULL, 1);
    return h;
  }

void SOGrid_Tree_free(SOGrid_Tree *t)
  { if (t != NULL)
    { dg_free_subtree(t->root);
      free(t);
    }
  }

#define SOGrid_Tree_FileFormat "2003-05-07"
  
void SOGrid_Tree_write(FILE *wr, SOGrid_Tree *t)
  { filefmt_write_header(wr, "SOGrid_Tree", SOGrid_Tree_FileFormat);
    fprintf(wr, "dim = %d\n", t->d);
    SOGrid_subtree_write(wr, t->root);
    filefmt_write_footer(wr, "SOGrid_Tree");
    fflush(wr);
  }
  
void SOGrid_subtree_write(FILE *wr, dg_tree_node_t *p)
  { if (p == NULL)
      { fputc('*', wr); }
    else
      { fputc('(', wr);
        SOGrid_subtree_write(wr, LOCH(*p)); 
        fputc(',', wr);
        SOGrid_subtree_write(wr, HICH(*p)); 
        fputc(')', wr);
      }
  }    
  
SOGrid_Tree *SOGrid_Tree_read(FILE *rd)
  { int dim;
    SOGrid_Tree *tree;
    filefmt_read_header(rd, "SOGrid_Tree", SOGrid_Tree_FileFormat);
    dim = nget_int(rd, "dim"); fget_eol(rd);
    affirm(dim >= 0, "bad dimension");
    tree = SOGrid_Tree_new(dim);
    fget_skip_formatting_chars(rd);
    tree->root = SOGrid_subtree_read(rd, NULL, 1);
    fget_skip_formatting_chars(rd);
    filefmt_read_footer(rd, "SOGrid_Tree");
    return tree;
  }
  
dg_tree_node_t *SOGrid_subtree_read(FILE *rd, dg_tree_node_t *p, dg_cell_index_t k)
  { int c = fget_char(rd);
  if (c == '*')
    { return NULL; }
  else if (c == '(')
    { dg_tree_node_t *this = dg_tree_node_new(p, k);
      // Read low subtree:
      LOCH(*this) = SOGrid_subtree_read(rd, this, 2*k);
      // Check for comma:
      fget_skip_formatting_chars(rd);
      c = fget_char(rd);
      affirm(c == ',', "comma expected");
      fget_skip_formatting_chars(rd);
      // Read high subtree:
      HICH(*this) = SOGrid_subtree_read(rd, this, 2*k + 1);
      // Check for close parenthesis:
      fget_skip_formatting_chars(rd);
      c = fget_char(rd);
      affirm(c == ')', " expected ')'");
      return this;
    }
  else
    { affirm(FALSE, "expected '*' or '('");
      return NULL;
    }
  }
