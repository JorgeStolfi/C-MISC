/* See {gch_solution.h} */
/* Last edited on 2014-07-23 19:50:25 by stolfilocal */

#define _GNU_SOURCE
#include <assert.h>
#include <stdlib.h>

#include <gem.h>

#include <gch_list.h>
#include <gch_problem_input.h>

#include <gch_solution.h>

void gch_solution_free(void *_solution)
  { gch_solution_t *solution = (gch_solution_t*)_solution;
    gch_problem_input_free(solution->pi);
    gch_list_free(solution->facets);
    free(solution);
  }

gem_ref_t gch_solution_duplicate(gch_solution_t *solution, int dim, gch_list_t *facets)
  { gem_ref_vec_t nodes = gem_ref_vec_new(1000);
    int i,j;
    int nn = 0;
    gem_traverse(solution->gem, dim, &nodes, &nn);
    for (i = 0; i < nn; i++)
      { assert(gem_get_label(nodes.e[i]) == i);
        gem_ref_t node = gem_node_new(dim);
        gem_set_data(node, gem_get_data(nodes.e[i]));
        for (j = 0; j < dim; j++)
          { gem_ref_t q = gem_step(nodes.e[i], j);
            int ind = gem_get_label(q);
            if (ind < i) { gem_splice(node, nodes.e[ind], j); }
          }
        nodes.e[i] = node;
      }
    gch_list_node_t *p = solution->facets->first;
    while(p != NULL)
      { gem_ref_t q = (gem_ref_t)(p->content);
        int ind = gem_get_label(q);
        gch_list_add(facets, nodes.e[ind]);
        p = p->next;
      }
    return nodes.e[gem_get_label(solution->gem)];
  }

int gch_solution_compare(void *a, void *b)
  {
    return gch_problem_input_compare(a, ((gch_solution_t*)b)->pi);
  }
