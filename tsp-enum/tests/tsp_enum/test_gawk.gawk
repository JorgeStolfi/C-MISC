#! /usr/bin/gawk -f
# Last edited on 2004-12-28 02:29:30 by stolfi

BEGIN {
  Pi = 3.141592653589793;         # An obscure constant sometimes used in statistics.
  for (i = 1; i <= 7; i++)
    { ti = sphere_slice_vol_fraction(i,0.5);
      ui = slice_pos_from_vol_fraction(i,ti,0.4,0.001);
      printf "%2d %24.16f %24.16f %24.16f\n", \
        i, sphere_vol(i), ti, ui \
        > "/dev/stderr";
    }
}

function slice_pos_from_vol_fraction(d,t,uini,utol,  ustep,umin,umax,u)
{
  # Finds {u} such that {sphere_slice_vol_fraction(d,u) = t}
  # Starts looking near {uini}.
  # Parameter {utol} is a tolerance for {u}.
  ustep = 3*utol;
  umin = uini;
  while (sphere_slice_vol_fraction(d,umin) > t) 
    { umin -= ustep ; if (umin < -1.0) { umin = -1.0; break; } }
  umax = uini;
  while (sphere_slice_vol_fraction(d,umax) < t) 
    { umax += ustep ; if (umax > 1.0) { umax = 1.0; break; } }
  while (1)
    { u = (umin + umax)/2;
      printf "  umin = %24.16f umax = %24.16f u = %24.16f", umin, umax, u > "/dev/stderr";
      if ((umax - umin) <= 2*utol) { break; }
      tu = sphere_slice_vol_fraction(d,u); 
      printf " tu = %24.16f\n", tu > "/dev/stderr";
      if (t < tu) 
        { umax = u; }
      else
        { umin = u; }
      u = (umin + umax)/2;
    }
  printf "\n" > "/dev/stderr";
  return u;
}

function sphere_slice_vol_fraction(d,u,  t)
{
  # Returns the fraction {T(d,u)} of the volume of the unit {d}-ball
  # that is contained in the slice between {x=-1} and {x=u},
  # for {u} in {[-1,+1]}.
  if (d < 0) { prog_error("bad d"); }
  if (d == 0)
    { return (u == -1.0 ? 0.25 : (u == 1 ? 0.75 : 0.5)); }
  else
    { t = atan2(u,sqrt(1-u*u)); # t = arcsin(u)
      return 0.5 + sphere_slice_vol_ang(d,t)/sphere_vol(d);
    }
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
        { C = exp((d-1)*log(C)); }
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

