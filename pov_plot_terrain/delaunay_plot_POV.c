/* See {delaunay_plot.h}. */
/* Last edited on 2009-01-06 04:22:27 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include <jsstring.h>
#include <quad.h>
#include <bool.h>

#include <delaunay.h>
#include <delaunay_plot_POV.h>


void delaunay_plot_POV_triangles(FILE *wr, quad_arc_t e,double height[],char* texture)
  {
    
    auto void draw_delaunay_triangle(quad_arc_t e);
    auto void draw_delaunay_triangle_pair(quad_arc_t e, void *closure);
    
    quad_enum(e, draw_delaunay_triangle_pair, NULL);
    
     void draw_delaunay_triangle_pair(quad_arc_t e, void *closure){
	draw_delaunay_triangle(e);
	draw_delaunay_triangle(quad_sym(e));
     }

    void draw_delaunay_triangle(quad_arc_t e)
      { delaunay_site_t *a = quad_org(e);
        delaunay_site_t *b = quad_dst(e);
	delaunay_site_t *c = quad_dst(quad_lnext(e));
	if( (a->index > b->index) || (a->index > c->index ) ){ return ; }
	double az,bz,cz;
	if( orient( a,b,c) <= 0 ){ 
		//Draw triangle on the reference plane
		az = bz = cz = 0;
	} 
	else{
		//draw traingle on surface
		az = height[a->index];
		bz = height[b->index];
		cz = height[c->index];	
	}
	
        
	fprintf(wr,"    triangle{ ");
	fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", a->p.c[0],a->p.c[1], az);
	fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", b->p.c[0],b->p.c[1], bz);
 	fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>",  c->p.c[0],c->p.c[1], cz);
	fprintf(wr,"  texture{ %s }",texture);
	fprintf(wr," }\n");
	return ;
      }

 }

void delaunay_plot_POV_skirt(FILE *wr, quad_arc_t e,double height[],char* texture)
  {
    
    auto bool_t triangle_is_CCW(quad_arc_t e);
    auto void draw_silouette_edge(quad_arc_t e);
    auto void draw_delaunay_triangle_pair(quad_arc_t e, void *closure);
    
    quad_enum(e, draw_delaunay_triangle_pair, NULL);
    
     void draw_delaunay_triangle_pair(quad_arc_t e, void *closure){
	bool_t L = triangle_is_CCW(e);
	bool_t R = triangle_is_CCW(quad_sym(e));
	if(L != R) draw_silouette_edge(e);
     }

    bool_t triangle_is_CCW(quad_arc_t e)
      { delaunay_site_t *a = quad_org(e);
        delaunay_site_t *b = quad_dst(e);
	delaunay_site_t *c = quad_dst(quad_lnext(e));
	return orient(a,b,c) > 0;
      }

      void draw_silouette_edge( quad_arc_t e){
	delaunay_site_t *a = quad_org(e);
        delaunay_site_t *b = quad_dst(e);	
	double az = height[a->index];
 	double bz = height[b->index];
	if (bz != 0){
		fprintf(wr,"    triangle{ ");
		fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", a->p.c[0],a->p.c[1], az);
		fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", b->p.c[0],b->p.c[1], bz);
		fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>",  b->p.c[0],b->p.c[1], 0.0);
		fprintf(wr,"  texture{ %s } ",texture);
		fprintf(wr," }\n");
	}
	if(az != 0 ){
		fprintf(wr,"    triangle{ ");
		fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", a->p.c[0],a->p.c[1], az);
		fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>,", b->p.c[0],b->p.c[1], 0.0);
		fprintf(wr,"  <%+9.6f,%+9.6f,%+9.6f>",  a->p.c[0],a->p.c[1], 0.0);
		fprintf(wr,"  texture{ %s } ",texture);
		fprintf(wr," }\n");
	}


      }
 }

