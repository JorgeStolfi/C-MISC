/* Last edited on 2014-07-22 02:36:04 by stolfilocal */
#ifndef gch_list_H
#define gch_list_H

#define _GNU_SOURCE

typedef struct gch_list_node_t 
  { void *content;
    struct gch_list_node_t *next, *prev;
  } gch_list_node_t;
  /* A node in a doubly-linked list. */

gch_list_node_t* gch_list_node_new(void *content);
void gch_list_node_free(gch_list_node_t *node);

typedef struct gch_list_t 
  { gch_list_node_t *first;
  } gch_list_t;
  /* Head of a doubly-linked list. */

typedef gch_list_node_t* gch_list_iterator_t;

gch_list_t* gch_list_new(void);
void gch_list_copy(gch_list_t *dest, gch_list_t *orig);
void gch_list_free(gch_list_t *list);
void gch_list_add(gch_list_t *list, void *content);
void gch_list_remove(gch_list_t *list, gch_list_iterator_t iter);
void gch_list_reset(gch_list_t *list, gch_list_iterator_t *iter);
void gch_list_advance(gch_list_iterator_t *iter);
int  gch_list_at_end(gch_list_iterator_t iter);
int  gch_list_empty(gch_list_t *list);
void* gch_list_content(gch_list_iterator_t iter);
void* gch_list_first_content(gch_list_t *list);
void* gch_list_first_extract(gch_list_t *list);

#endif
