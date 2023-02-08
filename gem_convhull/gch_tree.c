/* See {gch_tree.h} */
/* Last edited on 2014-07-22 02:43:56 by stolfilocal */

#include <stdlib.h>

#include <gch_tree.h>
#include <gch_util.h>

#define gch_tree_height(X) ((X) == NULL ? -1 : (X)->height)

gch_tree_t* gch_tree_new(void *content)
  { gch_tree_t *treegem_node_t = (gch_tree_t*)malloc(sizeof(gch_tree_t));
    treegem_node_t->content = content;
    treegem_node_t->height = 0;
    treegem_node_t->left = NULL;
    treegem_node_t->right = NULL;
    return treegem_node_t;
  }

void gch_tree_free(gch_tree_t *tree)
  { if (tree != NULL)
      { gch_tree_free(tree->left);
        gch_tree_free(tree->right);
        free(tree);
      }
  }

/* Balanced tree manipulation: */
gch_tree_t* gch_tree_simple_rotate_left(gch_tree_t *c);
gch_tree_t* gch_tree_simple_rotate_right(gch_tree_t *c);
gch_tree_t* gch_tree_internal_rotate_left(gch_tree_t *c);
gch_tree_t* gch_tree_internal_rotate_right(gch_tree_t *c);
gch_tree_t* gch_tree_balance(gch_tree_t *tree);

gch_tree_t* gch_tree_simple_rotate_left(gch_tree_t *c)
  { gch_tree_t *r = c->right;
    c->right = r->left;
    r->left = c;

    c->height = MAX(gch_tree_height(c->left),gch_tree_height(c->right))+1;
    r->height = MAX(gch_tree_height(r->left),gch_tree_height(r->right))+1;

    return r;
  }

gch_tree_t* gch_tree_simple_rotate_right(gch_tree_t *c)
  { gch_tree_t *l = c->left;
    c->left = l->right;
    l->right = c;

    c->height = MAX(gch_tree_height(c->left),gch_tree_height(c->right))+1;
    l->height = MAX(gch_tree_height(l->left),gch_tree_height(l->right))+1;

    return l;
  }

gch_tree_t* gch_tree_internal_rotate_left(gch_tree_t *c)
  { gch_tree_t *l = c;
    gch_tree_t *r = c->right;
    c = r->left;

    l->right = c->left;
    r->left = c->right;
    c->left = l;
    c->right = r;

    l->height = MAX(gch_tree_height(l->left),gch_tree_height(l->right))+1;
    r->height = MAX(gch_tree_height(r->left),gch_tree_height(r->right))+1;
    c->height = MAX(gch_tree_height(c->left),gch_tree_height(c->right))+1;

    return c;
  }

gch_tree_t* gch_tree_internal_rotate_right(gch_tree_t *c)
  { gch_tree_t *r = c;
    gch_tree_t *l = c->left;
    c = l->right;

    r->left = c->right;
    l->right = c->left;
    c->right = r;
    c->left = l;

    r->height = MAX(gch_tree_height(r->left),gch_tree_height(r->right))+1;
    l->height = MAX(gch_tree_height(l->left),gch_tree_height(l->right))+1;
    c->height = MAX(gch_tree_height(c->left),gch_tree_height(c->right))+1;

    return c;
  }

gch_tree_t* gch_tree_balance(gch_tree_t *tree)
  { int balance = gch_tree_height(tree->right) - gch_tree_height(tree->left);
    if (balance > 1)
      { if (gch_tree_height(tree->right->left) <= gch_tree_height(tree->right->right))
          { tree = gch_tree_simple_rotate_left(tree); }
        else
          { tree = gch_tree_internal_rotate_left(tree); }
      }
    else if (balance < -1)
      { if (gch_tree_height(tree->left->right) <= gch_tree_height(tree->left->left))
          { tree = gch_tree_simple_rotate_right(tree);  }
        else
          { tree = gch_tree_internal_rotate_right(tree); }
      }
    return tree;
  }

gch_tree_t* gch_tree_add(gch_tree_t *tree, void *content, int (*compare)(void*,void*))
  { if (tree == NULL)
      { return gch_tree_new(content); }
    else
      { if (compare(tree->content, content) > 0)
          { tree->left = gch_tree_add(tree->left, content, compare);
            tree->height = MAX(tree->height,tree->left->height+1);
          }
        else
          { tree->right = gch_tree_add(tree->right, content, compare);
            tree->height = MAX(tree->height,tree->right->height+1);
          }
        return gch_tree_balance(tree);
      }
  }

void* gch_tree_find(gch_tree_t *tree, void *signature, int (*compare)(void*,void*))
  { if (tree == NULL) { return NULL; }
    int comp = compare(tree->content,signature);
    if (comp < 0) { return gch_tree_find(tree->right,signature,compare); }
    if (comp > 0) { return gch_tree_find(tree->left,signature,compare); }
    return tree->content;
  }

void gch_tree_gem_traverse(gch_tree_t* tree, void (*action)(void*))
  { if (tree != NULL)
      { gch_tree_gem_traverse(tree->left,action);
        action(tree->content);
        gch_tree_gem_traverse(tree->right,action);
      }
  }

