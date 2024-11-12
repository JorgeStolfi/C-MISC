#! /usr/bin/gawk -f
# Last edited on 2013-10-26 00:01:00 by stolfilocal

BEGIN \
  {
    # A /spline pulse/ is a univariate polynomial spline {F} which has
    # knots at every integer and half-integer argument,
    # finite support, and is an even function (that is {F(x) = F(-x)}
    # for all {x}.
    # 
    # This program tries to determine the minimum support width {2*w} of a
    # spline pulse with each specified continuity order {c >= 0} and
    # degree {g} that is interpolative (that is, has {F(x) = 0} at every
    # nonzero integer {x}, and {F(0) = 1}) and can interpolate polynomials
    # of degree {g} exactly.

    cmin = 1;

    for (c = cmin; c <= 4; c++)
      { for (g = c+1; g <= 2*c+3; g++)
          { printf "c = %2d  g = %2d", c, g > "/dev/stderr";
            w = find_smallest_width(c, g, wmin);
            if (w == 1) { break; }
          }
        printf "\n" > "/dev/stderr";
      }
  }    


function find_smallest_width(c,g, wmin,w,ns,nc,nu,ne_c,ne_z,ne_i,ne_r,ne)
  {
    # Find the smallest pulse width {w} that has continuty {c} and degree {g}.
    wmin = (c == -1 ? 1 : (c == 0 ? 2 : 3));
    ns = 0; # Number of solutions for current {c,g}. 
    for (w = wmin; w < 3*(g+5); w++) 
      { # Parameters for positive half of grid.
        # Assumes the negative half is symmetric.

        # How many cells of width 1/2 on the positive half?
        nc = w;
        # How many unknowns (coefficients) in those cells?
        nu = nc*(g + 1);
        # How many continuity equations in the joints of those cells (except at {x = 0})?
        ne_c = nc*(c + 1);
        # How many continuity equations at {x = 0}?
        # All even derivatives are OK by symmetry;
        # the odd derivatives must be set to zero:
        ne_z = int((c + 1)/2);
        # How many interpolation equations?
        # We must have {F(x)=0} at every integer {x} except 0,
        # and {F(0)=1}.  The continuity equations already ensure that {F(x)=0}
        # at every joint {x >= nc/2}, so: 
        ne_i = int((nc + 1)/2);
        # How many equations to ensure reconstruction of 
        # degree-g polynomials? We already can get interpolation
        # at integer args, to get it for all args we need
        # {g} additional equations:
        ne_r = g;

        # How many equations total?
        ne = ne_c + ne_z + ne_i + ne_r;

        # Check solvability:
        if (nu >= ne)
          { # How many independent solutions?
            ns = nu - ne + 1;
            printf "  w = %2d  ns = %2d\n", w, ns > "/dev/stderr";
            break;
          }
        else
          { ns = 0; }
      }
    if (ns == 0) { printf " ** failed\n" > "/dev/stderr"; exit(1); }
  }
