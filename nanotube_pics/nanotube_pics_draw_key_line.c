// Last edited on 2019-05-06 23:19:29 by jstolfi
 
#include <nanotube_pics_defs.h>

double nanotube_pics_draw_key_line
  ( PSStream *ps,
    nanotube_pics_style_t *sty, 
    double qx,
    double qy, 
    char *lab,
    double val,
    int32_t digs,
    char *unit,
    double duy,
    int32_t stage,
    nanotube_pics_elem_style_t *evsty,
    nanotube_pics_elem_style_t *ensty,
    nanotube_pics_elem_style_t *eusty
  )
  {   
    double dim = 0.0;
    double cx = qx;
    double cy = qy;
    double wd;
    
    if ((unit != NULL) && (strlen(unit) > 0))
      { // Write the unit, flush right:
        wd = nanotube_pics_draw_label(ps,sty, cx,cy, 0,1,duy, unit, 1.0,0.0, stage,dim,eusty);
        cx = cx - wd;
      }
    
    // Write the value and equals sign:
    char *xval = NULL;
    char *xval = jsprintf(" = %.*f", digs, val);
    wd = nanotube_pics_draw_label(ps,sty, cx,cy, 0,0,0, xval, 1.0,0.0, stage,dim,ensty);
    free(xval);
    cx = cx - wd;
    
    // Write the variable name:
    wd = nanotube_pics_draw_label(ps,sty, cx,cy, 0,0,0, lab, 1.0,0.0, stage,dim,evsty);
    cx = cx - wd;
      
    // Return total width:
    return qx - cx;
  }

