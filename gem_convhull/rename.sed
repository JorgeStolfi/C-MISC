#! /bin/sed -f
# Last edited on 2014-07-22 02:15:32 by stolfilocal

# s:[	]: :g
# s:[)] *[{]:)\n  {:g
# s:for[(]:for (:g
# s:if[(]:if (:g
# s:[ ]*\([!]*[=]\)[ ]*: \1 :g
# s:[;] *:; :g
# s:[ ][ ]*$::g
s:[ ]*\([<>]\)[ ]*\([=]*\)[ ]*: \1\2 :g
s:[ ]*[-][ ]*[>][ ]*:->:g
s:[ ]*[=][ ]*[=][ ]*: == :g
# s:[ ]*\([-+*/|&]\)[ ]*[=][ ]*: \1= :g
/usage/s:[<][ ]*:<:g
/usage/s:[ ]*[>]:>:g

# s:[#]include *["]\(.*[.]h\)["] *$:#include <\1>:g
s:[#]include *[<] *\(.*[.]h\) *[>] *$:#include <\1>:g

# s:\b\(point_proj_hc\)\b:gch_\1:g
# s:\b\(pov_generator_hc\)\b:gch_\1_BUG:g
# s:\b\(pov_generator_gen\)\b:gch_pov_generator:g
# s:\b\(tree\)\b:gch_\1:g
# s:\b\(hyperplane\)\b:gch_\1:g
# s:\b\(point_generator_sphere\)\b:gch_\1:g
# s:\b\(basistransf\)\b:gch_\1:g
# s:\b\(util\)\b:gch_\1:g
# s:\b\(convexhull\)\b:gch_\1:g
# s:\b\(computedsolution\)\b:gch_\1:g
# s:\b\(giftwrap\)\b:gch_\1:g
# s:\b\(point_generator\)\b:gch_\1:g
# s:\b\(vector\)\b:gch_\1:g
# s:\b\(list\)\b:gch_\1:g
# s:\b\(vectorlist\)\b:gch_\1:g
# s:\b\(matrix\)\b:gch_\1:g
# s:\b\(geometry\)\b:gch_\1:g

# s:\bVectorList\b:gch_vector_list_t:g
# s:\bnewVectorList\b:gch_vector_list_read:g
# s:\bclearVectorList\b:gch_list_free:g
# s:\bVector\b:gch_vector_t:g
# s:\bnewVector\b:gch_vector_new:g
# s:\bsortIntArray\b:gch_int_array_sort:g
# 
# s:\bnewBasisTransf\b:gch_basistranf_new:g
# s:\bclearBasisTransf\b:gch_basistranf_free:g
# s:\btransfBasis\b:gch_basistranf_map:g
# s:\bBasisTransf\b:gch_basistranf_t:g
s:basistransf:affine_map:g
s:basistranf:affine_map:g
# 
# s:\bgetProjMatrix\b:gch_point_proj_matrix_make:g
# 
# # --> gch_problem_input
# s:\bnewProblemInput\b:gch_problem_input_new:g
# s:\bclearProblemInput\b:gch_problem_input_free:g
# s:\bcompareProblemInput\b:gch_problem_input_compare:g
# s:\bProblemInput\b:gch_problem_input_t:g
# 
# # --> gch_computed_solutions
# s:\bnewComputedSolutions\b:gch_computed_solutions_new:g
# s:\bclearComputedSolutions\b:gch_computed_solutions_free:g
# s:\bComputedSolutions\b:gch_computed_solutions_t:g
# s:\baddSolution\b:gch_computed_solutions_add:g
# s:\bcompareSolution\b:gch_solution_compare:g
# s:computedsolution[s]*:computed_solutions:g

# # --> gch_solution
# s:\bduplicateSolution\b:gch_solution_duplicate:g
# s:\bclearSolution\b:gch_solution_free:g
# s:\bfindSolution\b:gch_solution_find:g
# s:\bSolution\b:gch_solution_t:g
# 
# # --> gch_convexhull
# s:\bfindConvexHull\b:gch_convexhull_find:g
# s:\btrivialConvexHull\b:gch_convexhull_trivial:g
# s:\bfindConvexHull_rec\b:gch_convexhull_find_rec:g
# 
# # --> gch_geometry
# s:\bfindAdjacentFacet\b:gch_geometry_find_adj_facet:g
# s:\bfindFirstFacet\b:gch_geometry_find_first_facet:g
# s:\bgetAffineHullBasis\b:gch_geometry_affine_hull_basis_get:g
# 
# # --> gch_hyperplane
# s:\bnewHyperplane\b:gch_hyperplane_new:g
# s:\bclearHyperplane\b:gch_hyperplane_free:g
# s:\binitHyperplaneFromMatrix\b:gch_hyperplane_from_matrix:g
# s:\binitHyperplaneFromPoints\b:gch_hyperplane_from_points:g
# s:\btestHyperplanePoint\b:gch_hyperplane_point_test:g
# s:\btestHyperplanePointSign\b:gch_hyperplane_point_test_sign:g
# s:\bchangeHyperplaneSign\b:gch_hyperplane_neg:g
# s:\bHyperplane_print\b:gch_hyperplane_print:g
# s:\bHyperplane_isNull\b:gch_hyperplane_is_null:g
# s:\bHyperplane\b:gch_hyperplane_t:g
# 
# # --> gch_list
# s:\bnewList_node_t\b:gch_list_node_new:g
# s:\bclearList_node_t\b:gch_list_node_free:g
# s:\bnewList\b:gch_list_new:g
# s:\bcopyList\b:gch_list_copy:g
# s:\bclearList\b:gch_list_free:g
# s:\bList_add\b:gch_list_add:g
# s:\bList_remove\b:gch_list_remove:g
# s:\bList_reset\b:gch_list_reset:g
# s:\bList_advance\b:gch_list_advance:g
# s:\bList_atEnd\b:gch_list_at_end:g
# s:\bList_empty\b:gch_list_empty:g
# s:\bList_content\b:gch_list_content:g
# s:\bList_firstContent\b:gch_list_first_content:g
# s:\bList_extractFirst\b:gch_list_first_extract:g
# s:\bListIterator\b:gch_list_iterator_t:g
# s:\bList_node_t\b:gch_list_node_t:g
# s:\bList\b:gch_list_t:g
# 
# # --> gch_matrix
# s:\bnewMatrix\b:gch_matrix_new:g
# s:\bclearMatrix\b:gch_matrix_free:g
# s:\bsetRowVector\b:gch_matrix_set_row:g
# s:\bsubDet\b:gch_matrix_sub_det:g
# s:\bdet\b:gch_matrix_det:g
# s:\bMatrix_multVector\b:gch_matrix_map_col:g
# s:\bMatrix_print\b:gch_matrix_print:g
# s:\bMATRIX_ELEM\b:gch_matrix_elem:g
# s:\bMatrix\b:gch_matrix_t:g
# 
# # --> gch_tree
# s:\bsimpleLeftRotate\b:gch_tree_simple_rotate_left:g
# s:\bsimpleRightRotate\b:gch_tree_simple_rotate_right:g
# s:\binternalLeftRotate\b:gch_tree_internal_rotate_left:g
# s:\binternalRightRotate\b:gch_tree_internal_rotate_right:g
# s:\bTree_balance\b:gch_tree_balance:g
# s:\bclearTree\b:gch_tree_free:g
# s:\bTree_add\b:gch_tree_add:g
# s:\bTree_find\b:gch_tree_find:g
# s:\bTree_gem_traverse\b:gch_tree_gem_traverse:g
# s:\bTree\b:gch_tree_t:g
# 
# # --> gch_vector_list
# s:\bnewVectorList\b:gch_vector_list_new:g
# s:\bclearVectorList\b:gch_vector_list_free:g
# s:\bnewVectorSubList_Complete\b:gch_vector_list_sub_complete_new:g
# s:\bnewVectorSubList_CompleteProj\b:gch_vector_list_sub_complete_proj_new:g
# s:\bnewVectorSubList_Proj\b:gch_vector_list_sub_proj_new:g
# s:\bclearVectorSubList\b:gch_vector_list_sub_free:g
# s:\bgetVector_realIndex\b:gch_vector_list_real_index:g
# s:\bVectorSubList\b:gch_vector_list_sub_t:g
# s:\bVectorList\b:gch_vector_list_t:g
s:vectorlist:vector_list:g
# 
# # --> gch_vector
# s:\bnewVector\b:gch_vector_new:g
# s:\bclearVector\b:gch_vector_free:g
# s:\bVector_getCopy\b:gch_vector_copy:g
# s:\bVector_print\b:gch_vector_print:g
# s:\bVector\b:gch_vector_t:g




