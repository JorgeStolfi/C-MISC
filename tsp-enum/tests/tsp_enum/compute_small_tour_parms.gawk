#! /usr/bin/gawk -f
# Last edited on 2004-12-27 23:20:49 by stolfi

BEGIN {
  usage = ( "compute-small-tour-parms \\\n" \
    "  -v nv=NV [-v fixexpt=BOOL] < TPFILE > PARMS" \
  );
  abort = -1;
  
  # Reads a list of the {m} cheapest total tour weights, sorted,
  # and saves them in {Z[1..m]}. Adjusts the parameters
  # {Zref,Coef,Expt} so that the approximate formula 
  # {T(Z) = Coef*(Z - Zref)^Expt} is best satisfied.
  # Prints {Zref} and {Coef}.  If {fixexpt} is TRUE
  # sets {Expt = nv/2-1} instead of adjusting it.
  
  if (nv == "") { arg_error("must define nv"); }
  nread = 0;
  
  split("", Z);
  
  # Global constants for minimizer:
  Eps = (1.0e-4);                   # Square root of relative machine precision.
  Phi = (1.61803398874989484821);   # Golden ratio.
  Ihp = (0.61803398874989484821);   # 1/Phi.
  IhpSqr = (0.3819660112501051518); # 1/Phi^2.
  
  # Default/initial, min, and max {Expt}:
  minExpt = 0;
  maxExpt = nv;
  defaultExpt = nv/2.0;
  
  # Set to 1 to obtain debugging printouts:
  debug = 0;
}

(abort >= 0) { exit abort; }

// { Z[FNR] = $1; nread = FNR; next; }

END {
  if (abort >= 0) { exit abort; }
  ndata = nread;
  if (ndata == 0)
    { arg_error("not enough data"); }
  else if (ndata < 3)
    { # Not enough data for fitting, print a guess:
      printf "not enough data for proper fitting\n" > "/dev/stderr";
      Zref = 0.95*Z[1]; Coef = nv; Expt = defaultExpt;
    }
  else if (fixexpt)
    { 
      # Set fixed {Expt}, fit {Zref,Coef} to {Z[1]}--{Z[ndata]}:
      fit_zref_coef(defaultExpt);
    }
  else
    { 
      # Fit {Zref,Coef,Expt} to {Z[1]}--{Z[ndata]}:
      fit_expt_zref_coef();
    }
    
  # Output fitted params:
  printf "%24.16e %24.16e %24.16e\n", Zref, Coef, Expt;
  if (debug) 
    { printf "Zref = %24.16e Coef = %24.16e Expt = %24.16e\n", Zref, Coef, Expt \
        > "/dev/stderr";
    }
}

function fit_expt_zref_coef(  \
  xini,xa,xb,tol,dist \
)
{
  # Assumes {Z[1..ndata]} is set to the costs of the cheapest {ndata} 
  # solutions. Computes {Coef,Zref,Expt} that gives the best fit to the 
  # equation {i = Coef*(Z[i] - Zref)**Expt}
  
  xini = defaulExpt;
  xa = minExpt; xb = maxExpt;
  tol = 0.001;
  dist = xb-xa;
  brent_min(xini, tol, dist, xa, xb);
  
  # Set {Expt} and compute {Zref,Coef}, just to be sure:
  fit_zref_coef(x);
}

function brent_min( \
  xini,tol,dist,xa,xb, \
  u,fu,v,fv,w,fw, \
  d,e,r,Dfx,Dq,DDfx,DDq,fxOld,tol1,tol2,tol3 \
)
{
  # Looks for the minimum point of the globaly defined function 
  # {evalf}.  Requires a guess {xini}, a domain interval 
  # {xa,xb}, the desired accuracy of the answer {tol},
  # and an upper bound for the distamce between {xini} and the 
  # true minimum point. Sets the global variables {x,fx,dfx} 
  # to the solution.
  
  # v, fv: Third-smallest point of current triple
  # w, fw: Second-smallest point of current triple
  # u, fu: Splitting point

  # d, e, r:      Step sizes 
  # Dfx,Dq:       Estimated derivative at {x} is {Dfx/Dq}. 
  # DDfx,DDq:     Estimated second derivative at {x} is {DDfx/DDq}. 
  # fxOld:        Previous value of {fx} 
  
  tol3 = tol/3.0; 

  if ((xini < xa) || (xini > xb)) { prog_error("bad xini"); }
  
  x = xini; fx = evalf(x);
  v = x; fv = fx;
  w = x; fw = fx;

  e = 0.0;
  fxOld = fw;
  while (1)
  {
    if ((fw < fx) || (fw > fv)) { prog_error("bad fw"); }
    tol1 = Eps * abs(x) + tol3;
    tol2 = 2.000001 * tol1;
    if (debug) { printf "tol1 = %18.10f tol2 = %18.10f\n", tol1, tol2 > "/dev/stderr"; }

    if (debug) 
      { print_status(xa,xb,v,fv,x,fx,w,fw); }

    # check stopping criterion: 
    if (debug) { printf "x-xa =  %18.10f xb-x = %18.10f\n", x-xa,xb-x > "/dev/stderr"; }
      
    if ((x - xa <= tol2) && (xb - x <= tol2)) { return; }

    if (abs(e) > tol1)
      { # Fit parabola through {v}, {x}, {w}, sets {Dfx,Dq, DDfx,DDq}. 
        compute_derivatives(v, fv, w, fw, x, fx);
        r = e; e = d;
      }
    else
      { Dfx = 0.0; Dq = 0; DDfx = 0.0; DDq = 0.0;
        r = 0.0;
      }

    # Compute the displacement {d} from {x} to the next probe: 
    if (DDfx < 0.0) { Dfx = -Dfx; DDfx = -DDfx; }
    if ( \
       (abs(Dfx)*DDq < 0.5*abs(DDfx*r)*Dq) && \
       (Dfx*DDq < DDfx*(x - xa)*Dq) && \
       (Dfx*DDq > DDfx*(x - xb)*Dq) \
    ) 
      { # A parabolic-interpolation step: 
        if (debug) { printf "parabolic step\n" > "/dev/stderr"; }
        d = -(Dfx/DDfx)*(DDq/Dq);
      }
    else
      { # A golden-section interpolation: 
        if (debug) { printf "golden section\n" > "/dev/stderr"; }
        e = (xb - x > x - xa ? xb - x : xa - x);
        d = IhpSqr * e;
      }

    # The function must not be evaluated too close to {xa} or {xb}: 
    u = x + d;
    
    if (debug) { printf "u =  %18.10f\n", u > "/dev/stderr"; }
    
    if (u - xa < tol1)
      { u = xa + tol1; prte("A"); }
    else if (xb - u < tol1)
      { u = xb - tol1; prte("B"); }

    # The function must not be evaluated too close to {x}: 
    if (abs(u - x) >= tol1)
      { prte("C"); }
    else if (u > x)
      { u = x + tol1; prte("D"); }
    else
      { u = x - tol1; prte("E"); }  

    fu = evalf(u);

    if (debug) { printf "u =  %18.10f fu = %18.10f\n", u, fu > "/dev/stderr"; }
    
    # update  {xa,xb}: 
    if (fu >= fx)
      { if (u < x) { xa = u; } else { xb = u; } }
    if (fu <= fx)
      { if (u < x) { xb = x; } else{ xa = x; } }

    # Update {x,v,w}: 
    if (u < xa) { prog_error("u < xa"); }
    if (u > xb) { prog_error("u > xb"); }
    if (fu <= fx)
      { if ((v != w) && (w != x)) { v = w; fv = fw; }
        w = x; fw = fx;
        x = u; fx = fu;
      }
    else if (fu <= fw)
      { if (w == x) { prog_error("w == x"); } 
        v = w; fv = fw;
        w = u; fw = fu;
      }
    else if (fu <= fv)
      { if (v == x) { prog_error("v == x"); } 
        if (w == x)
          { w = u; fw = fu; }
        else
          { v = u; fv = fu; }
      }
    else if (v == w)
      { v = u; fv= fu; }
  }
}

function compute_derivatives ( \
  u,fu,v,fv,x,fx,
  uv, ux,xv,du,dv,q \
)
{
  # Computes estimates {Dfx/Dq} and {DDfx/DDq} for the first and second
  # derivatives of the goal function {F} at {x}, by fitting 
  # a quadratic through the three points {(u, fu)}, {(v, fv)},
  # and {(x, fx)}.  

  # The scaling factor {q} will be positive if all three arguments
  # {u}, {v}, {w} are sufficiently distinct, and zero otherwise.
  
  uv = v - u;
  ux = x - u;
  xv = v - x;
  du = xv * (fx - fu);
  dv = ux * (fv - fx);
  q = uv * ux * xv;
  
  Dfx = du*xv + dv*ux;
  DDfx = 2.0 * (dv - du);
  if (q < 0.0) { Dfx = -Dfx; DDfx = -DDfx; q = -q; }
  Dq = q;
  DDq = q;
}

function evalf(x)
{
  # Returns the total squared discrepancy assuming eponent {x}.
  # Warning: changes the global values of {Coef,Zref}.
  
  # Compute {Zref,Coef} for {Expt = x}: 
  fit_zref_coef(x);
  # Compute the goal function:
  return compute_discrepancy(Zref,Coef,x);
}

function fit_zref_coef(Ex,  i,Zi,Ti,Wi,SZZ,SZ1,S11,STZ,ST1,D,DZ,D1,KZ,K1)
{
  # Sets the global variable {Expt} to {Ex}, and
  # {Zref,Coef} to the best-fitting parameters for
  # eponent {Expt}.
  
  Expt = Ex;
  
  # Compute the normal system.
  # Basis: {Z} and {1}.
  # Scalar product { <f|g> = SUM{ f(Z[i])g(Z[i]) : i = 1..ndata }}.
  SZZ = 0; SZ1 = 0; S11 = 0;
  STZ = 0; ST1 = 0;
  for (i = 1; i <= ndata; i++)
    { Wi = 1.0;
      Zi = Z[i]; 
      Ti = exp(log(i)/Ex);
      SZZ += Wi*Zi*Zi; SZ1 += Wi*Zi; S11 += Wi;
      STZ += Wi*Ti*Zi; ST1 += Wi*Ti;
    }
  # Solve normal system {((SZZ, SZ1),(SZ1,S11)) * (CZ,C1) = (STZ,ST1)}: 
  D = SZZ*S11 - SZ1*SZ1;
  DZ = STZ*S11 - ST1*SZ1;
  D1 = SZZ*ST1 - SZ1*STZ;
  KZ = DZ/D;
  K1 = D1/D;
  if (debug) { printf "KZ = %24.16e D1 = %24.16e\n", KZ, K1 > "/dev/stderr"; }
  if (KZ < 0.0) { KZ = 0.001; }
  Coef = exp(log(KZ)*Ex);
  Zref = -K1/KZ;
  # Paranoid check:
  if (Zref <= 0.0) { Zref = 0.95*Z[1]; }
}

function compute_discrepancy(Zr,Cf,Ex,  S,i,Wi,Zest,Zobs,dZ)
{
  # Computes the total square discrepancy between {(i/Cf)**(1/Ex)}
  # and {Z[i]-Zr}, weighted by {i}.
  S = 0;
  for (i = 1; i <= ndata; i++)
    { Wi = 1.0;
      Zest = exp(log(i/Cf)/Ex) + Zr;
      Zobs = Z[i]; 
      dZ = Zest-Zobs
      S += Wi*dZ*dZ;
    }
  return S;
}

function print_status(xa,xb,v,fv,x,fx,w,fw)
{
  printf "\n" > "/dev/stderr";
  printf "xa = %18.10f\n", xa > "/dev/stderr";
  printf "v =  %18.10f fv = %18.10f\n", v, fv > "/dev/stderr";
  printf "x =  %18.10f fx = %18.10f\n", x, fx > "/dev/stderr";
  printf "w =  %18.10f fw = %18.10f\n", w, fw > "/dev/stderr";
  printf "xb = %18.10f\n", xb > "/dev/stderr";
  printf "\n" > "/dev/stderr";
}

function abs(x)
{
  return (x < 0 ? -x : x);
}

function prte(tag)
{
  if (debug) { printf "[%s]\n", tag > "/dev/stderr"; } 
}

function arg_error(msg)
{
  printf "** %s\n", msg > "/dev/stderr";
  printf "usage: %s\n", usage > "/dev/stderr";
  abort = 1;
  exit abort;
}

      
      
