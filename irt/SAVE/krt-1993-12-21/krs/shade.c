#include <math.h>
#include "shade.h"
#include "raymath.h"

void shade(
    t_3d *pos, t_3d *ray, t_3d *nrm,
    t_object *obj,
    t_color *color
)
{
    int lnum;
    double k, dis, bright, spec, diffuse;
    t_surface *surf;
    t_3d refl, ltray;
    t_color newcol;

    /* calculate reflected ray */
    k = -2.0 * dotp( ray, nrm );
    refl.x = k * nrm->x + ray->x;
    refl.y = k * nrm->y + ray->y;
    refl.z = k * nrm->z + ray->z;

    /* ambient light contribuition */
    surf = &surface[obj->surfnum];
    color->r = surf->ar;
    color->g = surf->ag;
    color->b = surf->ab;

    for( lnum=0; lnum < nlight; lnum++ ){

	/* get ray to light */
	lightray( lnum, pos, &ltray );
	diffuse = dotp( nrm, &ltray );
	if( diffuse > 0.0 ){

	    /* object faces light, add diffuse */
	    bright = brightness( obj->id, lnum, pos, &ltray );
	    diffuse *= bright;
	    color->r += surf->dr * diffuse;
	    color->g += surf->dg * diffuse;
	    color->b += surf->db * diffuse;

	    spec = dotp( &refl, &ltray );
	    if( spec > 0.0 ){

		/* highlight is here, add specular */
		spec = bright * pow( spec, surf->coef );
		color->r += surf->sr * spec;
		color->g += surf->sg * spec;
		color->b += surf->sb * spec;
	    }
	}
    }

    /* reflection */
    k = surf->refl;
    if( ( k > 0.0 ) && ( level < maxlevel ) ){

	level++;
	dis = intersect( obj->id, pos, &refl, &newcol );
	if( dis > 0 ){
	    color->r += newcol.r * k;
	    color->g += newcol.g * k;
	    color->b += newcol.b * k;
	}
	else{
	    color->r += background.r * k;
	    color->g += background.g * k;
	    color->b += background.b * k;
	}
	level--;
    }

    /* transparency */
    k = surf->transp;
    if( k > 0.0 ){

	color->r *= ( 1-k );
	color->g *= ( 1-k );
	color->b *= ( 1-k );
	dis = intersect( obj->id, pos, ray, &newcol );
	if( dis > 0 ){
	    color->r += newcol.r * k;
	    color->g += newcol.g * k;
	    color->b += newcol.b * k;
	}
	else{
	    color->r += background.r * k;
	    color->g += background.g * k;
	    color->b += background.b * k;
	}
    }
}

#include "light.h"
#include "raymath.h"

double s_litdis;

void lightray(
    int lnum,
    t_3d *objpos, t_3d *lray
)
{
    lray->x = light[lnum].x - objpos->x;
    lray->y = light[lnum].y - objpos->y;
    lray->z = light[lnum].z - objpos->z;
    s_litdis = normalize( lray );
}

double brightness(
    int source, int lnum,
    t_3d *pos, t_3d *ray
)
{
    int iobj;
    double s;

    for( iobj=0; iobj < nobject; iobj++ ){

	if( iobj != source ){ /* don't try source */
	    s = intsolid ( pos, ray, object[iobj].solid );
	    if( ( s > 0.0 ) && ( s <= s_litdis ) )
	        return( 0.0 );  /* object in sahdow */
	}
    }
    /* object not in shadow */
    return( light[lnum].bright );
}
