/* Last edited on 2014-07-22 02:39:03 by stolfilocal */
#ifndef gch_tree_H
#define gch_tree_H
/* Balanced search trees. */

#define _GNU_SOURCE

typedef struct gch_tree_t
  { void *content;
    int height;
    struct gch_tree_t *left, *right;
  } gch_tree_t;

gch_tree_t* gch_tree_new(void *content);

void gch_tree_free(gch_tree_t *tree);

gch_tree_t* gch_tree_add(gch_tree_t *tree, void *content, int (*compare)(void*,void*));

void* gch_tree_find(gch_tree_t *tree, void *signature, int (*compare)(void*,void*));

void gch_tree_gem_traverse(gch_tree_t* tree, void (*action)(void*));

#endif
