
/* Last edited on 2020-01-01 09:26:23 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <values.h>
#include <string.h>
#include <assert.h>

#include <sheet_cut.h>

#include <vec.h>
#include <argparser.h>
#include <r2.h>
#include <bool.h>
#include <fget.h>
#include <affirm.h>

#include <csc_group.h>

/* INTERNAL PROTOTYPES  */

sheet_cut_node_t* csc_find_and_delete_plate(char *tag, sheet_cut_node_vec_t *plv);
  /* Assumes that {plv.e[0..plv.ne-1]} is a list of {sheet_cut_node_t]
    pointers, some of them maybe {NULL}.  Finds the node with the given {tag}
    and returns it, after setting that entry of {plv} to {NULL}. 
    
    Fails if there is no such node in {plv}, or there are two nodes with the 
    same tag. */

void csc_group_read_one_item(FILE *rd, int32_t *cmdP, char **tagP);
  /* Reads the next item from {rd}, in the format {csc_GROUPING_INFO} below.
    
    If the next item in {rd} is a command code, ('[', ']', '@', etc.),
    stores that character into {*cmdP}, and {NULL} into {*tagP}.  
    
    If the next item is a plate tag, stores that tag into {*tagP}, and 0 in {*cmdP}.
    
    if it runs into end-of-file, returns {EOF = -1} in {*cmdP} and {NULL} 
    in {*tagP}.  */
    
sheet_cut_node_t *csc_group_nodes
  ( sheet_cut_node_vec_t *itv, 
    double cut_wd
  );
  /* Packs the items (plates or blocks) {itv.e[0..itv.ne-1]} as
    children of a new block node.  The items are packed horizontally,
    with their bottom sides at ordinate 0, starting at abscissa 0.
    The procedure leaves a gap of width {cut_wd} between them.

    The items {itv.e[0..itv.ne-1]} must all have {.pos} at {(0,0)}. 
    The field {itv.e[i].pos} will be set by the
    procedure to the assigned position of item {itv.e[i]} 
    relative to the low corner of the block.  
    The block itself will have {.pos =(0,0)}.
    
    When the procedure starts or exits, the block {st.blk} may
    or may not be empty, and the {st}'s canvas may be empty, full,
    or only partially occupied. */


/* IMPLEMENTATIONS */

sheet_cut_node_vec_t csc_group_read_items
  ( FILE *rd, 
    double cut_wd,
    sheet_cut_node_vec_t *plv,
    int32_t level
  )
  { 
    /* Debugging stuff: */
    bool_t debug = FALSE;
    char *dbpref = NULL;
    if (debug) { asprintf(&dbpref, "%*scsc_read_and_pack_items[%d]:", 2*level, "", level); }
    
    /* Vector of top-blocks created: */
    sheet_cut_node_vec_t topv = sheet_cut_node_vec_new(20); 
    int32_t ntop = 0; /* Current top blocks are {topv.e[0..ntop-1]}. */

    /* Loop on input items (piece descriptions or commands) until EOF or ']': */
    bool_t flip_next = FALSE; /* TRUE to flip the next plate or block. */
    while (TRUE)
      { int32_t cmd;
        char *tag;  
        csc_group_read_one_item(rd, &cmd, &tag);
        if (debug) 
          { fprintf(stderr, "%s cmd = %d", dbpref, cmd);
            if (cmd >= ' ') { fprintf(stderr, " = '%c'", cmd); }
            fprintf(stderr, " tag = %s\n", (tag == NULL ? "NULL" : tag));
          }
        if ((cmd == 0) || (cmd == '['))
          { /* New plate or group {pl}. */
            sheet_cut_node_t *it = NULL;
            if (cmd == 0)
              { /* A single plate. */
                /* Find tag in plate list, and remove it from there: */
                it = csc_find_and_delete_plate(tag, plv);
                assert(it != NULL);
              }
            else
              { /* A sub-group. Recursively read the items up to the ']' and pack them into a block: */
                sheet_cut_node_vec_t subv = csc_group_read_items(rd, cut_wd, plv, level+1);
                /* Require the ']': */
                fget_skip_formatting_chars(rd);
                int32_t ch = fgetc(rd);
                if (ch != ']') 
                  { fprintf(stderr, "**error: missing ']'\n");
                    assert(FALSE);
                  }
                /* Pack the items as a block: */
                it = csc_group_nodes(&subv, cut_wd);
              }
            if (flip_next) { sheet_cut_flip_node(it); }
            sheet_cut_node_vec_expand(&topv, ntop);
            topv.e[ntop] = it;
            ntop++;
            flip_next = FALSE;
          }
        else if (cmd == EOF) 
          { /* End of file: */ 
            break;
          }
        else if (cmd == ']')
          { /* End of a group started by caller; leave the ']' to be parsed by it: */
            ungetc(cmd, rd);
            break;
          }
        else if (cmd == '@')
          { /* Flip (or unflip) the next item: */
            flip_next = TRUE;
          }
        else 
          { assert(FALSE); }
      }
    sheet_cut_node_vec_trim(&topv, ntop);
    return topv;
  }
  
void csc_group_read_one_item(FILE *rd, int32_t *cmdP, char **tagP)
  { 
    bool_t debug = FALSE;
    (*cmdP) = 0;  /* "No command". */
    (*tagP) = NULL; /* "No plate". */
    while (TRUE)
      { /* Read another line, which may be a comment: */
      fget_skip_formatting_chars(rd);
      int32_t ch = fgetc(rd);
      if (debug) 
        { fprintf(stderr, "csc_read_item: ch = %d", ch);
          if (ch >= ' ') { fprintf(stderr, " = '%c'", ch); }
          fputc('\n', stderr);
        }
      if (ch == EOF)
        { /* No more input, return: */
          (*cmdP) = EOF;
          return;
        }
      else if (ch == '#')
        { /* Ignore rest of line and try again: */
          fget_skip_to_eol(rd); 
        }
      else if ((ch == '[') || (ch == ']') || (ch =='@'))
        { /* Command code, return it: */
          (*cmdP) = ch;
          return;
        }
      else
        { /* Must be a plate tag: */
          ungetc(ch, rd);
          char *tag = fget_string(rd);
          assert(tag != NULL);
          if (debug) { fprintf(stderr, "  tag = %s", tag); }
          /* Return the plate just read: */
          (*tagP) = tag;
          return;
        }
    }
  }

sheet_cut_node_t *csc_group_nodes
  ( sheet_cut_node_vec_t *itv,
    double cut_wd
  )
  {
    int32_t nit = itv->ne; /* Number of items to group. */
    demand(nit > 0, "group cannot be empty");
    sheet_cut_node_t* blk = NULL;
    double cur_cx = 0; /* Abscissa of next child. */
    for (int32_t kit = 0; kit < nit; kit++)
      { sheet_cut_node_t *itk = itv->e[kit];
        assert(itk != NULL);
        assert((itk->pos.c[0] == 0.0) &&(itk->pos.c[1] == 0.0));
        if (blk != NULL) { cur_cx += cut_wd; }
        itk->pos = (r2_t){{ cur_cx, 0.0 }};
        blk = sheet_cut_add_child(blk, itk);
        cur_cx += itk->size.c[0];
      }
    assert(blk != NULL);
    return blk;
  } 

sheet_cut_node_t* csc_find_and_delete_plate(char *tag, sheet_cut_node_vec_t *plv)
  {
    sheet_cut_node_t *plr = NULL;
    for (int32_t k = 0; k < plv->ne; k++)
      { sheet_cut_node_t *plk = plv->e[k];
        if ((plk != NULL) && (strcmp(plk->tag, tag) == 0))
          { if (plr != NULL) 
              { fprintf(stderr, "** two plates with tag = \"%s\"\n", tag); exit(1); }
            plr = plk; plv->e[k] = NULL;
          }
       }
    if (plr == NULL)
      { fprintf(stderr, "** no plate with tag = \"%s\", or used twice\n", tag); exit(1); }
    return plr;
  }
