// Last edited on 2019-05-13 02:58:39 by jstolfi

#include <nanotube_pics_defs.h>

void nanotube_pics_fig_strip_choose_a1
  ( nanotube_pics_style_t *sty,
    double a0x,
    double a0y,
    double ux,
    double uy,
    double vx,
    double vy,
    double krx,
    double kry,
    double wx,
    double wy,
    double *a1xP,
    double *a1yP
  )
  {
    auto double eval_a1(double a1x, double a1y);
      /* Evaluates the candidate placement {a1x,a1y} for the left end
        of {w}. Lower values mean better choices. Returns a number
        between 0 and 1000 if all constraints are satisfied; between
        1000 and 2000 if the only violation is that the strip overlaps
        the key; and 2000 or higher if other constraints are
        violated. */
    
    double dx = sty->dx;
    double dy = sty->dy;
    double xytot = hypot(sty->xtot, sty->ytot); // biggest dist inside piture
    
    // Minimum clearances between arrow ends and things, to fit the labels "A1","A2":
    double fnsize = sty->labsty.fontsize;
    double mrgax = 1.2*fnsize*mm; // Min hor clearance
    double mrgay = 1.0*fnsize*mm; // Min ver clearance
    
    double mrgs = 0.5*dx; // Min hor clearance between strip and things.

    // Index ranges worth considering:
    int32_t k1xmin = (int32_t)ceil((mrgax - sty->orgx)/dx);
    int32_t k1xmax = (int32_t)floor((sty->xtot - mrgax - wx - sty->orgx)/dx);
    affirm(k1xmin <= k1xmax, "not wide enough");
    
    int32_t k1ymin = (int32_t)ceil((mrgay - sty->orgy)/dy);
    int32_t k1ymax = (int32_t)floor((sty->ytot - mrgay - wy - sty->orgy)/dy);
    affirm(k1ymin <= k1ymax, "not tall enough");
    
    // Try all positions and save the best:
    (*a1xP) = NAN; (*a1yP) = NAN;
    double evmin = +INF; // Best evaluation so far.
    for (int32_t k1y = k1ymin; k1y <= k1ymax; k1y++)
      for (int32_t k1x = k1xmin; k1x <= k1xmax; k1x++)
        { double a1x = sty->orgx + (k1x + 0.5*(k1y % 2))*dx;
          double a1y = sty->orgy + k1y*dy;
          
          double ev = eval_a1(a1x, a1y);
          if (ev < evmin) { (*a1xP) = a1x; (*a1yP) = a1y; evmin = ev; }
        }
    fprintf(stderr, "best a1 = (%.1f,%.1f) qual = %12.5f\n", ((*a1xP)-sty->orgx)/dx, ((*a1yP)-sty->orgy)/dy, evmin); 
    return;

    // Internal implementations:

    double eval_a1(double a1x, double a1y)
      { 
        if (wx < 0)
          { // Evaluation is the same as if reversed:
            a1x = a1x + wx; wx = -wx;
            a1y = a1y + wy; wy = -wy;
          }
          
        double wL = hypot(wx, wy); // Length of arrow.
        
        // Right end of arrow:
        double a2x = a1x + wx;
        double a2y = a1y + wy;

        // The quality is how close the center of {w} is to the center of the picture:
        double evx = (a1x + a2x)/2 - sty->xtot/2;
        double evy = (a1y + a2y)/2 - sty->ytot/2;
        double qual = 100*(evx*evx + evy*evy)/(xytot*xytot);
        
        // Arrow and labels lie withing the margins:
        assert (a1x <= a2x);
        bool_t soft_ok = TRUE; // Satisfies all soft constraints.
        bool_t hard_ok = TRUE; // Satisfies all hard constraints.
        
        // Arrow and labels must fit inside margins:
        hard_ok &= ((a1x >= mrgax) && (a2x <= sty->xtot - mrgax));
        hard_ok &= ((fmin(a1y,a2y) >= mrgay) && (fmax(a1y,a2y) <= sty->ytot - mrgay));
        
        // Tip of {u} must be {mrgs} away from strip:
        double upos = ((a0x + ux - a1x)*wx + (a0y + uy - a1y))/(wL*wL); 
        hard_ok &= (upos <= 0 - mrgs/wL) && (upos >= 1 + mrgs/wL);

        // Tip of {v} must be {mrgs} away from strip:
        double vpos = ((a0x + vx - a1x)*wx + (a0y + vy - a1y))/(wL*wL);
        hard_ok &= (vpos <= 0 - mrgs/wL) && (vpos >= 1 + mrgs/wL);

        // Arrow and labels must not overlap {u,v} frame:
        if (wy >= 0)
          { // Assume that {u v} frame is at bot left:
            hard_ok &= ((a1x >= a0x + ux + mrgax) || (a1y >= a0y + vy + mrgay));
          }
        else
          { // Assume that {u v} frame is at top left:
            hard_ok &= ((a1x >= a0x + ux + mrgax) || (a1y <= a0y - mrgay));
          }
        
        // Arrow and labels must not overlap the key:
        if (wy >= 0)
          { // Assume that key is at top left:
            hard_ok &= ((a2x <= krx - mrgax) || (a2y <= kry - mrgay));
          }
        else
          { // Assume that key is at bot left:
            hard_ok &= ((a2x <= krx - mrgax) || (a2y >= kry + mrgay));
          }
        
        // Corner of key must be {mrgs} away from strip (soft):
        double kpos = ((krx - a1x)*wx + (kry - a1y))/(wL*wL); 
        soft_ok &= (kpos <= 0 - mrgs/wL) && (kpos >= 1 + mrgs/wL);

        // Conclusion:
        if (! hard_ok)
          { return 2000 + qual; }
        else if (! soft_ok)
          { return 1000 + qual; }
        else
          { return qual; }
      }
      
  }
