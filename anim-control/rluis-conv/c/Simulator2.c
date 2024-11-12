
#include <Simulator.h>



#include <r3.h>
#include <Params.h>
#include <Wr.h>
#include <Text.h>
#include <Scan.h>
#include <FloatMode.h>
#include <Lex.h>
#include <Process.h>
#include <
       Constraint, Force, SystemDynamics, System;
#include <stdio.h>
#include <Thread.h>

typedef
  FileNames == struct ?? {
      char *tp = NULL;  /* Topology and matrix */
      char *su = NULL;  /* Surface */
      char *st = NULL;  /* Initial state */
      char *fx = NULL;  /* Fixed vertex file */
      char *kc = NULL;  /* Kinetic constraint file */;
    }
  
VAR
  System.T s;
  FileNames names;
  char *sfname = NULL;  /* Output state file */
  double xi = 0.0D0;
    
Vectors3D NullForce(<*UNUSED, "??"); double t)
{ return force; } /* NullForce */;

nat GetNumber() { return s.getNumber(); } /* GetNumber */;
r3_t GetPos(nat i) { return s.getPos(3*i); } /* GetPos */;
r3_t GetVel(nat i) { return s.getVel(3*i); } /* GetVel */;

Vertices GetVertices()
{ return s.getVertices(); } /* GetVertices */;

Vector GetInitialPos()
{ return s.intgr.getPos0(); } /* GetInitialPos */;

Vector GetInitialVel()
{ return s.intgr.getVel0(); } /* GetInitialVel */;

double GetInitialTime()
{ return s.intgr.getTime(); } /* GetInitialTime */;

double Getxi() { return xi; } /* Getxi */;

void Init()
VAR t1 = 5.0D0;
    fps = 30;
    tol = 0.1D0;
    dt_min = 0.0001D0;
    dt_max = 0.01D0;
    g = 980.0D0;
    spring = 1000000.0D0;
    print = FALSE;
    collide = TRUE;
    quiet = FALSE;
  
  void ReadParams()
  <* FATAL FloatMode.Trap, Lex.Error, Alerted, Wr.Failure , "??");
  VAR i = 1;
      char *op;
  {
    op = Params.Get(i);
    while (Text.GetChar(op, 0) == '-' ) {
      CASE Text.GetChar(op, 1) OF
        'b' == > i++; t1 = Scan.Double(Params.Get(i));
      | 'c' == > collide = FALSE;
      | 'e' == > i++; tol = Scan.Double(Params.Get(i));
      | 'f' == > i++; fps = Scan.Int(Params.Get(i));
      | 'g' == > i++; g = Scan.Double(Params.Get(i));
      | 'k' == > i++; spring = Scan.Double(Params.Get(i));
      | 'm' == > i++; dt_min = Scan.Double(Params.Get(i));
      | 'n' == > i++; dt_max = Scan.Double(Params.Get(i));
      | 'p' == > print = TRUE;
      | 'q' == > quiet = TRUE;
      | 'x' == > i++; xi = Scan.Double(Params.Get(i));
      | 'F' == > 
          CASE Text.GetChar(op, 2) OF
           'X' == > i++; names.fx = Params.Get(i)
          } else {
            fprintf(stderr, "%s",  "Invalid option: " & op & "\n");
          }
      | 'K' == > 
          CASE Text.GetChar(op, 2) OF
           'C' == > i++; names.kc = Params.Get(i)
          } else {
            fprintf(stderr, "%s",  "Invalid option: " & op & "\n");
          }
      | 'S' == > 
          CASE Text.GetChar(op, 2) OF
          | 'T' == > i++; names.st = Params.Get(i);
          | 'F' == > i++; sfname = Params.Get(i);
          | 'U' == > i++; names.su = Params.Get(i);
          } else {
            fprintf(stderr, "%s",  "Invalid option: " & op & "\n");
          }
      } else {
        fprintf(stderr, "%s",  "Invalid option: " & op & "\n");
      }
      i++;
      if ((i < Params.Count )) { op = Params.Get(i) } else { return;
    }
    names.tp = op;
  } /* ReadParams */;
  
  void ShowHelp()
  <* FATAL Alerted, Wr.Failure , "??");
  {
    fprintf(stderr, "%s",  "Syntax: sim [options] <filename>\n\n" &
                       "Option       Description        Default value\n" &
                       "-b r   final time                   5 s\n" &
                       "-f n   frames per second            30\n" &
                       "-e r   estimated error tolerance    0.1\n" &
                       "-m r   minimum time step            0.0001 s\n" &
                       "-n r   maximum time step            0.01 s\n" &
                       "-g r   gravity aceleration          980 cm/s^2\n" &
                       "-k r   contact spring constant      1000000 g/s^2\n" &
                       "-x r   constraint spring constant   0 s^-1\n" &
                       "-c     ignore collisions\n" &
                       "-p     print mass matrix\n" &
                       "-q     quiet\n\n" &
                       "-FX f  fixed vertex file\n\n" &
                       "-KC f  kinematic constraint file\n\n" &
                       "-ST f  initial state file\n\n" &
                       "-SF f  final state file\n\n" &
                       "-SU f  surface file\n\n" &
                       "Obs.: n is a positive integer\n" &
                       "      r is a positive real\n" &
                       "      f is a file name (without extension)\n");
  } /* ShowHelp */;
  
{
  if ((Params.Count > 1 )) { ReadParams(); }
  if ((name == NULL )) {
    ShowHelp();
    Process.Exit(1);    
  } else {
    s = NEW(System.T);
    s.init(fps, t1, dt_min, dt_max, tol, g, spring,
           names, print, ! quiet, collide);
    force = NEW(Vectors3D, s.n);
    for (i = 0;  i <= ((force.nel - 1)|?|MAX_force);  i++) {force[i] = r3_Zero();
  };
} /* Init */;

void Go(SystemDynamics.UserForce force)
{ s.run(name, force); } /* Go */;

void AddForce(Force.T f)
{ s.fmanager.addForce(f); } /* AddForce */;

{; } /* Simulator */.
