/* See {gch_list.h} */
/* Last edited on 2014-07-22 02:38:11 by stolfilocal */

#define _GNU_SOURCE
#include <stdlib.h>

#include <gch_list.h>

gch_list_node_t* gch_list_node_new(void *content)
  { gch_list_node_t *node = (gch_list_node_t*)malloc(sizeof(gch_list_node_t));
    node->content = content;
    node->next = NULL;
    node->prev = NULL;
    return node;
  }

void gch_list_node_free(gch_list_node_t *node)
  { free(node);
  }

gch_list_t* gch_list_new()
  { gch_list_t *list = (gch_list_t*)malloc(sizeof(gch_list_t));
    list->first = NULL;
    return list;
  }

void gch_list_copy(gch_list_t *dest, gch_list_t *orig)
  { gch_list_node_t* p = orig->first;
    while(p != NULL)
      { gch_list_add(dest,p->content);
        p = p->next;
      }
  }

void gch_list_free(gch_list_t *list)
  { gch_list_node_t* p = list->first;
    gch_list_node_t* q;
    while(p != NULL)
      { q = p->next;
        gch_list_node_free(p);
        p = q;
      }
    free(list);
  }

void gch_list_reset(gch_list_t *list, gch_list_iterator_t *iter)
  {
    *iter = list->first;
  }

void gch_list_advance(gch_list_iterator_t *iter)
  {
    *iter = (*iter)->next;
  }

int  gch_list_at_end(gch_list_iterator_t iter)
  {
    return (iter == NULL);
  }

int  gch_list_empty(gch_list_t *list)
  {
    return (list->first == NULL);
  }

void* gch_list_content(gch_list_iterator_t iter)
  {
    return iter->content;
  }

void* gch_list_first_content(gch_list_t *list)
  {
    return list->first->content;
  }

void* gch_list_first_extract(gch_list_t *list)
  { void *content = list->first->content;
    gch_list_remove(list,list->first);
    return content;
  }

void gch_list_add(gch_list_t *list, void *content)
  { gch_list_node_t* p = gch_list_node_new(content);
    if (list->first != NULL) { list->first->prev = p; }
    p->next = list->first;
    list->first = p;
  }

void gch_list_remove(gch_list_t *list, gch_list_iterator_t iter)
  { if (iter != NULL)
      { if (iter->prev != NULL) { iter->prev->next = iter->next; }
        if (iter->next != NULL) { iter->next->prev = iter->prev; }
        if (iter == list->first) { list->first = iter->next; }
        gch_list_node_free(iter);
      }
  }
