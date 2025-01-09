// Last edited on 2019-05-13 02:28:09 by jstolfi

#include <nanotube_pics_defs.h>

nanotube_pics_style_t *nanotube_pics_def_style(double xtot, double ytot)
  {
    nanotube_pics_style_t *sty = notnull(malloc(sizeof(nanotube_pics_style_t)), "no mem");

    sty->xtot = xtot; // Canvas width (mm).
    sty->ytot = ytot; // Canvas height (mm).
    
    sty->dy = 12*mm;
    sty->dx = sty->dy*2/sqrt(3);
    
    sty->bdlen = hypot(sty->dx/2, sty->dy/3);
    
    sty->b0x = 0.0;
    sty->b0y = 0.0;
    
    sty->orgx = sty->b0x + sty->dx/4;
    sty->orgy = sty->b0y + sty->dy/6;

    // General style parameters:
    sty->bgcol = (frgb_t){{ 0.930f, 0.930f, 0.930f }}; //  Background color.

    // Shadow displacements:
    sty->shcol =  (frgb_t){{ 0.667f, 0.667f, 0.667f }};  //  Color of shadows.
    sty->shx = -2*px; // Shadow displacement X.
    sty->shy = -3*px; // Shadow displacement Y.
    sty->rsh = 1.5*px;  // Extra width of shadow, on each side.
 
    // White halos:
    sty->whcol = sty->bgcol;  //  Color of halos.

    // Parameters for Latin variables in key:
    sty->kyvsty.fcol = (frgb_t){{ 0.000f, 0.300f, 0.200f }};
    sty->kyvsty.dcol = (frgb_t){{ 0.000f, 0.300f, 0.200f }};
    sty->kyvsty.lw  = 0.0;    // Chars are only filled,not stroked.  
    sty->kyvsty.rdot = 0.0;
    sty->kyvsty.arlen = 0.0;
    sty->kyvsty.arwid = 0.0;
    sty->kyvsty.rwh = 0.0;
    sty->kyvsty.font  = "ArialItalic"; 
    sty->kyvsty.fontsize  = 14;          

    // Parameters for Greek variables and symbols in key:
    sty->kyssty.fcol = (frgb_t){{ 0.000f, 0.300f, 0.200f }};
    sty->kyssty.dcol = (frgb_t){{ 0.000f, 0.300f, 0.200f }};
    sty->kyssty.lw  = 0.0;    // Chars are only filled,not stroked.  
    sty->kyssty.rdot = 0.0;
    sty->kyssty.arlen = 0.0;
    sty->kyssty.arwid = 0.0;
    sty->kyssty.rwh = 0.0;
    sty->kyssty.font  = "Symbol"; 
    sty->kyssty.fontsize  = 14;          

    // Parameters for numbers and equal signs in key:
    sty->kynsty.fcol = (frgb_t){{ 0.000f, 0.300f, 0.200f }};
    sty->kynsty.dcol = (frgb_t){{ 0.000f, 0.300f, 0.200f }};
    sty->kynsty.lw  = 0.0;    // Chars are only filled,not stroked.  
    sty->kynsty.rdot = 0.0;
    sty->kynsty.arlen = 0.0;
    sty->kynsty.arwid = 0.0;
    sty->kynsty.rwh = 0.0;
    sty->kynsty.font  = "Arial"; 
    sty->kynsty.fontsize  = 14;          

    // Parameters for general labels:
    sty->labsty.fcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->labsty.dcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->labsty.lw  = 0.0;    // Labels are only filled,not stroked.  
    sty->labsty.rdot = 0.0;
    sty->labsty.arlen = 0.0;
    sty->labsty.arwid = 0.0;
    sty->labsty.rwh = 3.0*px;
    sty->labsty.font  = "ArialItalic"; 
    sty->labsty.fontsize  = 22;          

    // Parameters for circumference vector and dots:
    sty->rfvsty.fcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->rfvsty.dcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->rfvsty.lw  = 3.0*px;         
    sty->rfvsty.rdot = 3*sty->rfvsty.lw;
    sty->rfvsty.arlen = 9*sty->rfvsty.lw;
    sty->rfvsty.arwid = 4*sty->rfvsty.lw;
    sty->rfvsty.rwh = 3.0*px;
    sty->rfvsty.font  = "ArialItalic"; 
    sty->rfvsty.fontsize  = 22;          

    // Parameters for strip and sector boundary lines:
    sty->edgsty.fcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->edgsty.dcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->edgsty.lw  = 2.5*px;         
    sty->edgsty.rdot = 3*sty->edgsty.lw;
    sty->edgsty.arlen = 9*sty->edgsty.lw;
    sty->edgsty.arwid = 4*sty->edgsty.lw;
    sty->edgsty.rwh = 3.0*px;
    sty->edgsty.font  = "ArialItalic"; 
    sty->edgsty.fontsize  = 22;          

    // Parameters for basis vector and dots:
    sty->bassty.fcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->bassty.dcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->bassty.lw  = 2.5*px;         
    sty->bassty.rdot = 3*sty->bassty.lw;
    sty->bassty.arlen = 7*sty->bassty.lw;
    sty->bassty.arwid = 3*sty->bassty.lw;
    sty->bassty.rwh = 3.0*px;
    sty->bassty.font  = "ArialItalic"; 
    sty->bassty.fontsize  = 22;          

    // Parameters for non-isomorphic vectors and dots in master image:
    sty->vatsty.dcol = (frgb_t){{ 0.930f, 0.133f, 0.000f }};
    sty->vatsty.fcol = (frgb_t){{ 0.930f, 0.133f, 0.000f }};
    sty->vatsty.lw  = 2.0*px;         
    sty->vatsty.rdot = 8*sty->vatsty.lw;
    sty->vatsty.arlen = 12*sty->vatsty.lw;
    sty->vatsty.arwid = 4*sty->vatsty.lw;
    sty->vatsty.rwh = 1.5*px;
    sty->vatsty.font  = "Arial"; 
    sty->vatsty.fontsize  = 12;          
    
    // Parameters enantiomeric vectors and dots:
    sty->vbtsty.dcol = (frgb_t){{ 0.233f, 0.000f, 0.930f }};
    sty->vbtsty.fcol = (frgb_t){{ 0.233f, 0.000f, 0.930f }};
    sty->vbtsty.lw  = 2.0*px;         
    sty->vbtsty.rdot = 8*sty->vbtsty.lw;
    sty->vbtsty.arlen = 12*sty->vbtsty.lw;
    sty->vbtsty.arwid = 4*sty->vbtsty.lw;
    sty->vbtsty.rwh = 2*px;
    sty->vbtsty.font  = "Arial"; 
    sty->vbtsty.fontsize  = 12;          
    
    // Parameters for type pair labels:
    sty->vntsty.dcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->vntsty.fcol = (frgb_t){{ 0.000f, 0.000f, 0.500f }};
    sty->vntsty.lw  = 0.0;         // Don't stroke chars.
    sty->vntsty.rdot = 0.0;
    sty->vntsty.arlen = 0.0;
    sty->vntsty.arwid = 0.0;
    sty->vntsty.rwh = 5*px;
    sty->vntsty.font  = "Arial"; 
    sty->vntsty.fontsize  = 14;          
 
    // Parameters for graphene and nanotube lattice:
    sty->bdsty.dcol = (frgb_t){{ 0.067f, 0.267f, 0.067f }};
    sty->bdsty.fcol = (frgb_t){{ 0.067f, 0.267f, 0.067f }};
    sty->bdsty.lw  = 2.5*px;    
    sty->bdsty.rdot = 0.0;
    sty->bdsty.arlen = 12*sty->bdsty.lw;
    sty->bdsty.arwid = 4*sty->bdsty.lw;
    sty->bdsty.rwh = 4.0*px;
    sty->bdsty.font  = "ArialItalic"; 
    sty->bdsty.fontsize  = 22;          

    for (uint32_t k = 0;  k < 2; k++)
      { 
        sty->atsty[k].lw  = 2.5*px;         
        sty->atsty[k].rdot = 30*px;
        sty->atsty[k].arlen = 0.0;
        sty->atsty[k].arwid = 0.0;
        sty->atsty[k].rwh = 4.0*px;
        sty->atsty[k].font  = "ArialItalic"; 
        sty->atsty[k].fontsize  = 22; 
      }
    sty->atsty[0].dcol = (frgb_t){{ 0.133f, 0.500f, 0.000f }};
    sty->atsty[0].fcol = (frgb_t){{ 0.333f, 1.000f, 0.000f }};

    sty->atsty[1].dcol = (frgb_t){{ 0.067f, 0.400f, 0.000f }};
    sty->atsty[1].fcol = (frgb_t){{ 0.200f, 0.800f, 0.000f }};
    
    return sty;
  }
