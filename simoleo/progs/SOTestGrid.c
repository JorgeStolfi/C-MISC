/* SOTestGrid.c */
/* Last edited on 2004-06-19 18:23:28 by stolfi */

#include <SOGrid.h>

#include <dg_grid.h>

#include <stdio.h>

int main (int argc, char **argv)
{
  interval_t box[2];
  dg_dim_t d = 2;
  dg_cell_index_t index;
  dg_rank_t r = 4;
  
  printf(" Cell Index: \n");
  scanf("%ld", (dg_cell_index_t*)&index);

  dg_cell_box_root_relative(index, r, d, box);

  printf(" Cell %ld with Rank %d and Dim %d \n", index, r, d);
 
  printf(" min_0 %6.4f max_0 %f \n", LO(box[0]), LO(box[0]));
  printf(" min_1 %6.4f max_1 %f \n", LO(box[1]), LO(box[1]));

  return 0;
}
