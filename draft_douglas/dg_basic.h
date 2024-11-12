#ifndef dg_basic_H
#define dg_basic_H

/* Basic definitions for multi-dimensional dyadic grids */
/* Last edited on 2014-05-15 22:49:23 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <psp_basic.h>
#include <psp_grid.h>

/* 
   SPACE DIMENSIONS AND COORDINATE AXES */

#define dg_dim_MAX (4)
#define dg_axis_MAX (dg_dim_MAX-1)
  /* Just to be safe. High grid dimensions {d} are dangerous because
    too many things have size {c^d} for some {c>1}.  In particular, 
    a {box_face_index_t} ranges from {0..3^d-1}. */

typedef psp_dim_t dg_dim_t; 
  /* Dimension of a box, space, etc.. */

typedef psp_axis_t dg_axis_t; 
  /* Identifies a coordinate axis of the enclosing space, from 0. */

typedef psp_axis_set_t dg_axis_set_t; 
  /* A packed subset of the coordinate axes. */

typedef psp_axis_index_t dg_axis_index_t; 
  /* Identifies a coordinate axis amonga set: 0 is the lowest. */

#define dg_axis_NONE (psp_axis_NONE)
  /* A {dg_axis_t} value that means "no axis". */

#define dg_axis_NOT_FOUND (psp_axis_NOT_FOUND)
  /* A {dg_axis_t} value that means "no such axis" */
  
dg_axis_set_t dg_axis_set_complement(dg_dim_t d, dg_axis_set_t A);
  /* The complement of {A} relative to the set {0..d-1}. */

/*
  CONTINUITY ORDER AND POLYNOMIAL DEGREE */
  
typedef psp_cont_t dg_cont_t;
  /* Continuity order of a function. The value {-1} means discontinuous.. */
  
typedef psp_degree_t dg_degree_t;
  /* Degree of a polynomial, polynomial spline, etc. along some axis. */

#endif
