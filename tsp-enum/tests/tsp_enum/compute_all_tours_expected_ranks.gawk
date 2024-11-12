#! /usr/bin/gawk -f
# Last edited on 2005-06-17 07:47:40 by stolfi

BEGIN {
  usage = ( "compute-all-tour-parms \\\n" \
    "  -v nv=NV [-v method=NUM] < TPFILE > PARMS" \
  );
  abort = -1;
  
  # Reads a list of the weights {Z[1..m]} of {m} random tours, sorted.
  # Adjusts the parameters {ZC,ZR,Dim} so that the approximate formula
  # {t(Z) = Sph(Dim, (Z - ZC)/ZR)} is best satisfied, where {t(Z)} is
  # the fraction of tours with cost at most {Z}, and {Sph(d,u)} is the
  # fraction of the volume of {d}-dimensional sphere of unit radius
  # between {x=-1} and {x=u}. Then writes out a list of pairs {Z[i],
  # t(Z[i])}
  # 
  # If {method=0} sets {Dim = nv-3}, {ZC} to the mean of {Z[1..m]}, 
  # and {ZR} based on the extremal values.
  # 
  # If {method=1} sets {Dim = nv-3}, then fits {ZC} and {ZR}
  # by least squares (treating {t} as the independent variable).
  # 
  # If {method=2} ajusts {Dim} too by nonlinear minimization.
  
  if (nv == "") { arg_error("must define nv"); }
  if (method == "") { arg_error("must define method"); }
  
  split("", Z); # Tour values, sorted in increasing order.
  
  # Global constants for minimizer:
  Eps = 1.0e-4;                   # Square root of relative machine precision.
  Phi = 1.61803398874989484821;   # Golden ratio.
  Ihp = 0.61803398874989484821;   # 1/Phi.
  IhpSqr = 0.3819660112501051518; # 1/Phi^2.
  Pi = 3.141592653589793;         # An obscure constant sometimes used in statistics.
  
  # Default/initial {Dim}, minimum {Dim}, maximum {Dim}
  defaultDim = 2*nv-1;
  minDim = defaultDim - 3; if (minDim < 1) { minDim = 1; }
  maxDim = defaultDim + 3;
  
  # Number of data points read:
  M = 0;
  
  # Set to 1 to obtain debugging printouts:
  debug = 1;
}

(abort >= 0) { exit abort; }

/^ *$/ { next; }

// { 
  M++;
  Z[M] = $1;
  if ((M > 1) && (Z[M] < Z[M-1]))
    { data_error(("data out of order")); }
  nread = FNR;
  next;
}

END {
  if (abort >= 0) { exit abort; }
  M = nread;
  printf "read %d tour values\n", M > "/dev/stderr";

  # Fit coefficients to data:
  if (method == 0)
    { # Quick fitting based on mean and extremal values:
      if (M < 2) { arg_error("not enough data"); }
      quick_fit_zc_zr(defaultDim);
    }
  else if (method == 1)
    { # Set fixed {Dim}, fit {ZC,ZR} to {Z[1]}--{Z[M]}:
      if (M < 2) { arg_error("not enough data"); }
      fit_zc_zr(defaultDim);
    }
  else
    { # Fit {ZC,ZR,Dim} to {Z[1]}--{Z[M]}:
      if (M < 3) { arg_error("not enough data"); }
      fit_dim_zc_zr();
    }
    
  # Display fitted params:
  if (debug) 
    { printf "ZC = %24.16e ZR = %24.16e Dim = %24.16e\n", ZC, ZR, Dim \
        > "/dev/stderr";
    }
    
  # Output pairs {Z[i],t(Z[i])}:
  for (i = 1; i <= M; i++)
    { ui = (Z[i] - ZC)/ZR;
      if (ui < -1.0) { ui = -1.0; }
      if (ui > 1.0) { ui = 1.0; }
      if (Dim < 0) { prog_error("bad Dim"); }
      tZi = sphere_slice_vol_fraction(Dim, ui);
      printf "%24.16e %24.16e\n", Z[i], tZi;
    }
}

function sphere_slice_vol_fraction(d,u,  t,f)
{
  # Returns the fraction {T(d,u)} of the volume of the unit {d}-ball
  # that is contained in the slice between {x=-1} and {x=u},
  # for {u} in {[-1,+1]}.
  # if (debug) { printf "    vol_fraction(%.4f,%.10f)", d, u > "/dev/stderr"; } 
  if (d < 0) { prog_error("bad d"); }
  if (d == 0)
    { f = (u == -1.0 ? 0.25 : (u == 1 ? 0.75 : 0.5)); }
  else
    { t = atan2(u,sqrt(1-u*u)); # t = arcsin(u)
      f = 0.5 + sphere_slice_vol_ang(d,t)/sphere_vol(d);
    }
  # if (debug) { printf " = %.10f\n", f > "/dev/stderr"; }
  return f;
}

function sphere_slice_vol_ang(d,t,   C,A,B)
{
  # Returns the volume of the slice of the unit {d}-ball
  # between {x=0} and {x=sin(t)},
  # for {t} in {[-Pi/2,Pi/2]}, namely 
  # {F(d,t) = V(d-1)*integral((cos(z))^d, z=0..t)}.
  if (d <= 0) { prog_error("bad d"); }
  if (d == 1)
    { return sin(t); }
  else if (d == 2)
    { return sin(t)*cos(t) + t; }
  else
    { C = cos(t); 
      if (C <= 0.0)
        { C = 0; }
      else if (C >= 1.0)
        { C = 1.0; }
      else
        { C = (d-1)*log(C); C = (C < -300 ? 0 : exp(C)); }
      A = sphere_vol(d-1)*sin(t)*C;
      B = 2*Pi*sphere_slice_vol_ang(d-2,t);
      return (A + B)/d;
    }
}

function sphere_vol(d)
{
  # Returns the volume of the unit {d}-ball, namely
  # {V(d) = Pi^(d/2)/((d/2)!)}
  if (d < 0) { prog_error("bad d"); }
  if (d == 0)
    { return 1; }
  else if (d == 1)
    { return 2; }
  else
    { return Pi/(d/2) * sphere_vol(d-2); }
}

function fit_zc_zr(Dm,  SUU,SU1,S11,SZU,SZ1,i,Utol,Ui,Ti,Wi,D,DU,D1,KU,K1)
{
  # Sets the global variable {Dim} to {Dm}, and
  # {ZC,ZR} to the best-fitting parameters for
  # that dimension.
  
  Dim = Dm;
  
  # Tolerance for {Ui} computation:
  Utol = (0.5/Dim)/M;
  
  # Compute the normal system, fitting {KU*U(i) + K1} to {Z[i]},
  # where {U(i)} is {u} such that {T(Dim,u) = (i-0.5)/M}.
  # Basis: {U(i)} and {1}
  # Scalar product { <f|g> = SUM{ f(i)g(i) : i = 1..M }}.
  SUU = 0; SU1 = 0; S11 = 0;
  SZU = 0; SZ1 = 0;
  Ui = -1.0;
  for (i = 1; i <= M; i++)
    { Wi = 1.0;
      Ti = (i-0.5)/M;
      # Compute {Ui} using previous value as initial guess:
      Ui = slice_pos_from_vol_fraction(Dim,Ti,Ui,Utol);
      Zi = Z[i]; 
      SUU += Wi*Ui*Ui; SU1 += Wi*Ui; S11 += Wi;
      SZU += Wi*Zi*Ui; SZ1 += Wi*Zi;
    }
  # Solve normal system {((SUU, SU1),(SU1,S11)) * (KU,K1) = (SZU,SZ1)}: 
  D = SUU*S11 - SU1*SU1;
  DU = SZU*S11 - SZ1*SU1;
  D1 = SUU*SZ1 - SU1*SZU;
  KU = DU/D;
  K1 = D1/D;
  if (debug) { printf "KU = %24.16e K1 = %24.16e\n", KU, K1 > "/dev/stderr"; }
  # Now compute {ZC,ZR}:
  ZC = K1;
  ZR = KU;
  # Paranoid checks:
  if (ZR <= 0) { ZR = 0.001; }
  if (ZR <= ZC - Z[1]) { ZR = 1.01*(ZC - Z[1]); }
  if (ZR <= Z[M] - ZC) { ZR = 1.01*(Z[M] - ZC); }
}

function slice_pos_from_vol_fraction(d,t,uini,utol,  ustep,umin,umax,u,tu)
{
  # Finds {u} such that {sphere_slice_vol_fraction(d,u) = t}
  # Starts looking near {uini}.
  # Parameter {utol} is a tolerance for {u}.
  ustep = 3*utol;
  umin = uini;
  while (1)
    { if (umin < -1.0) { umin = -1.0; break; }
      if (sphere_slice_vol_fraction(d,umin) <= t) { break; }
      umin -= ustep;
    }
  umax = uini;
  while (1)
    { if (umax > 1.0) { umax = 1.0; break; }
      if (sphere_slice_vol_fraction(d,umax) >= t) { break; }
      umax += ustep;
    }
  while (1)
    { u = (umin + umax)/2;
      # printf "  umin = %24.16f umax = %24.16f u = %24.16f", umin, umax, u > "/dev/stderr";
      if ((umax - umin) <= 2*utol) { break; }
      tu = sphere_slice_vol_fraction(d,u); 
      # printf " tu = %24.16f\n", tu > "/dev/stderr";
      if (t < tu) 
        { umax = u; }
      else
        { umin = u; }
      u = (umin + umax)/2;
    }
  # printf "\n" > "/dev/stderr";
  return u;
}

function fit_dim_zc_zr(  \
  d,S,dBest,SBest \
)
{
  # Assumes {Z[1..M]} is set to the costs of {M} 
  # solutions. Computes {ZC,ZR,Dim} that gives the best fit to the 
  # equation {i = ZC*(Z[i] - ZR)**Dim}
  
  dBest = minDim;
  SBest = evalf(dBest);
  for (d = minDim+1; d <= maxDim; d++)
    { S = evalf(d);
      if (S < SBest) { dBest = d; SBest = S; }
    }
  if (debug) { printf "  best Dim = %d\n", dBest > "/dev/stderr"; }
  # Set {Dim} and compute {ZC,ZR}, just to be sure:
  fit_zc_zr(dBest);
}

function evalf(x,  S)
{
  # Returns the total squared discrepancy assuming eponent {x}.
  # Warning: changes the global values of {ZC,ZR}.
  
  # Compute {ZC,ZR} for {Dim = x}: 
  printf "evalf(%.10f)", x > "/dev/stderr";
  fit_zc_zr(x);
  # Compute the goal function:
  S = compute_discrepancy(ZC,ZR,x);
  printf " = %.10f\n", S > "/dev/stderr";
  return S;
  
}

function compute_discrepancy(Zc,Zr,Dm,  SWZ,SW,i,Wi,Ti,Ui,Zest,Zobs,dZ,Utol)
{
  # Computes the total square discrepancy between {(i/Zr)**(1/Dm)}
  # and {Z[i]-Zc}, weighted by {i}.

  SWZ = 0;
  SW = 0;
  Ui = -1.0;
  Utol = (0.5/Dm)/M;
  for (i = 1; i <= M; i++)
    { Wi = 1.0;
      Ti = (i-0.5)/M;
      Ui = slice_pos_from_vol_fraction(Dm,Ti,Ui,Utol);
      Zest = Zc + Ui*Zr;
      Zobs = Z[i]; 
      dZ = Zest-Zobs
      SWZ += Wi*dZ*dZ;
      SW += Wi;
    }
  return SWZ/SW;
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

function data_error(msg)
{
  printf "%s:%s: ** %s\n", FILENAME, FNR, msg > "/dev/stderr";
  abort = 1;
  exit abort;
}

function prog_error(msg)
{
  printf "** %s\n", msg > "/dev/stderr";
  abort = 1;
  exit abort;
}

      
      
