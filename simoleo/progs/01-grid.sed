# Last edited on 2005-06-05 20:51:35 by stolfi

# s/SOGrid_rank/dg_cell_rank/g
# s/SOGrid_min_index/dg_cell_min_index/g
# s/SOGrid_max_index/dg_cell_max_index/g
# s/SOGrid_longest_axis/dg_split_axis/g
# s/SOGrid_cell_coords/dg_cell_box_root_relative/g
# s/SOGrid_Rank/dg_Rank/g
# s/SOGrid_Index/dg_CellIndex/g
# s/SOGrid_Height/dg_Rank/g
# s/SOGrid_Direction/dg_BinaryDir/g
# s/SOGrid_Dim/dg_Dim/g
# s/SOGrid_Axis/dg_Axis/g 
# s/SOGRID_MAX_DIM/DG_MAX_GRID_DIM/g

# min[k] -> box[k].end[0]
# max[k] -> box[k].end[1]
# SOGrid_cell_coords(k,r,d,min,max) -> dg_cell_box_root_relative(d,k,box)

# s/<dgbasic.h>/<dg_grid.h>/g
# s/<dgtree.h>/<dg_tree.h>/g

