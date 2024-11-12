MODULE SystemDynamics;

IMPORT Wr, LR3, LR3x3;
IMPORT SystemTopology, Integrator;
IMPORT 
  SparseSymetricMatrix AS SSM, 
  KineticConstraints AS KC,
  ViscoelasticForces AS VF, 
  ExternalForces AS EF,
  ContactForces AS CF,
  CollisionDetector AS CD;
  
FROM State IMPORT Coords;
FROM Thread IMPORT Alerted;
FROM Stdio IMPORT stderr;

REVEAL T = Public BRANDED OBJECT
    nCoords: CARDINAL;
  OVERRIDES
    getMassMatrix := GetMassMatrix;
    computeAcclerations := ComputeAccelerations
  END;

TYPE
  NodeMasses = ARRAY OF LONGREAL;

PROCEDURE GetMassMatrix(
    s: T; 
    mmRd: Rd.T;   (* Mass matrix file, or NIL if to be computed. *)
    mmWr: Wr.T;   (* Output factored mass matrix (may be NIL). *)
    verbose: BOOL := FALSE; 
  ) =
  <* FATAL Wr.Failure, Alerted *>
  BEGIN
    s.nCoords := 3 * s.nNodes;
    WITH
      nCoords = s.nCoords,
      fp = s.fingerprint()+BoolFingerprint(s.fixed)
    DO
      (* Mass matrix: *)
      WITH
        
      DO
        IF mmRd # NIL THEN
          IF verbose THEN Message("Reading factored mass matrix...") END;
          SSM.Read(mmRd, s.M, fp);
          IF verbose THEN Message("\n") END;
        ELSE
          IF verbose THEN
            Message("Computing factored mass matrix...");
            ResetCPUTime();
          END;
          s.M := ComputeFactoredMassMatrix(s, s.cf.fixed^);
          IF verbose THEN
            WriteCPUTime(stderr);
            Message("\n")
          END;
        END;
        IF mmWr # NIL THEN
          IF verbose THEN Message("Writing factored mass matrix...") END;
          SSM.Write(mmWr, s.M, fp);
          IF verbose THEN Message("\n") END;
        END
      END;

      (* Allocate force vector: *)
      s.force := NEW(REF Coords, nCoords);
      
    END;
  END Init;
  
PROCEDURE ComputeAcceleration(
    s: T;
    READONLY p: Position;    (* Generalized coordinates of system *)
    READONLY v: Velocity;    (* Time derivative of position. *)
    VAR f: Force;            (* IN: External forces; OUT: total forces. *)
    VAR J: ConstraintMatrix; (* Constraint gradients. *)
    VAR psi: Acceleration;   (* Constraint-relative acceleration. *)
    VAR D: DirectionMatrix;  (* Direction vectors OF reaction forces. *)
    VAR a: Acceleration;     (* OUT: Second derivatives of position. *)
    VAR h: Force;            (* OUT: Magnitude of reaction forces. *)
  ) =
  BEGIN
    <* ASSERT s.M # NIL *>
    AddInternalForces(s, p, v, f);
    (* Clear out forces on fixed coordinates: *)
    FOR k := 0 TO LAST(f) DO IF s.kc.fixed[k] THEN f[k] := 0.0d0 END END;
    (* Compute reaction forces for other constraints: *)
    AddReactionForces(s.M, f, J, psi, D, h);
    FOR k := 0 TO LAST(f) DO <* ASSERT f[k]*0.0d0 < 0.0d0 *> END;
    (* Compute accelerations: *)
    SSM.Solve(s.M, f, a);
    FOR k := 0 TO LAST(a) DO 
      <* ASSERT a[k]*0.0d0 < 0.0d0 *>
      IF s.kc.fixed[k] THEN <* ASSERT a[k] = 0.0d0 *> END
    END
  END ComputeAcceleration;

PROCEDURE ComputeFixedCells(
    s: SystemTopology.T; 
    RAEDONLY coordFixed: BOOLS;
  ) REF BOOLS =
  (*
    Returns a vector "cf" such that "cf[c]=TRUE" iff all corners of cell "c"
    are completely fixed. *)
    
  BEGIN
    WITH
      nCells = s.nCells,
      rcf = NEW(REF BOOLS, nCells), cf =  rcf^
    DO
      FOR c := 0 TO nCells-1 DO
        cf[c] := TRUE;
        WITH cell = s.cell[c] DO
          FOR j := 0 TO 3 DO
            WITH n = cell.node[j] DO
              FOR k := 3*n TO 3*n+2 DO 
                IF NOT s.kc.fixed[k] THEN cf[c] := FALSE END
              END
            END
          END
        END
      END;
      RETURN rcf
    END
  END ComputeFixedCells;

PROCEDURE ComputeFactoredMassMatrix(
    s: SystemTopology.T;
    READONLY coordFixed: BOOLS;
  ): SSM.T =
  <* FATAL Wr.Failure, Alerted *>
  BEGIN
    WITH
      M = AssembleMassMatrix(s)
    DO
      SSM.Remove(M, coordFixed);
      SSM.Factor(M);
      RETURN M
    END
  END ComputeFactoredMassMatrix;

PROCEDURE AssembleMassMatrix(s: SystemTopology.T): SSM.T =
  BEGIN
    WITH
      nNodes = s.nNodes,
      nCells = s.nCells,
      nCoords = s.nCoords,
      M = SSM.New(nCoords);
      row0 = NEW(REF Coords, nCoords),
      row1 = NEW(REF Coords, nCoords),
      row2 = NEW(REF Coords, nCoords)
    DO
      FOR n := 0 TO nNodes-1 DO
        FOR j := 0 TO nCoords-1 DO
          row0[j] := 0.0D0; 
          row1[j] := 0.0D0; 
          row2[j] := 0.0D0;
        END;
        WITH node = s.node[n], i = 3*n DO
          FOR k := 0 TO node.nCells DO
            WITH 
              c = node.cell[k],
              cell = s.cell[c]
            DO
              FOR v := 0 TO 3 DO
                WITH m = cell.node[v], j = 3*m DO
                  row0[j+0] := row0[j] + cell.mass/20.0D0;
                  row1[j+1] := row0[j];
                  row2[j+2] := row0[j];
                END
              END;
            END
          END;
          row0[i+0] := 2.0D0*row0[i];
          row1[i+1] := row0[i];
          row2[i+2] := row0[i];
          SSM.PutRow(M, i+0, row0);
          SSM.PutRow(M, i+1, row1);
          SSM.PutRow(M, i+2, row2);
        END
      END;
      RETURN M
    END
  END AssembleMassMatrix;

PROCEDURE ComputeNodeMasses(s: SystemTopology.T; gravity: LONGREAL): REF NodeMasses =
  BEGIN
    WITH
      nNodes = s.nNodes,
      nCells = s.nCells,
      rwt = NEW(REF Coords, nNodes),
      wt = rwt^
    DO
      FOR n := 0 TO nNodes-1 DO wt[n] := 0.0d0 END;
      IF gravity # 0.0d0 THEN
        FOR c := 0 TO nCells-1 DO
          WITH 
            cell = s.cell[c]
          DO
            FOR v := 0 TO 3 DO
              WITH n = cell.node[v] DO
                wt[n] := wt[n] + cell.mass/4.0D0;
              END
            END;
          END
        END;
        FOR n := 0 TO nNodes-1 DO wt[n] := gravity*wt[n] END;
      END;
      RETURN rwt
    END
  END ComputeNodeMasses;

PROCEDURE ComputeInternalForces(
    s: T;
    READONLY pos, vel: Coords;
  ) =
  VAR f := s.f;
      g := s.intgr;
      fx := s.fixed;
  BEGIN
    WITH
      f = s.force^,
      nodeFx = s.nodeFixed^
    DO
      FOR c := 0 TO s.nCells-1 DO
        IF NOT cellFixed[c] THEN
          WITH
            cell = s.cell[c],
            node = cell.node,
            k0 = 3*node[0], 
            k1 = 3*node[1], 
            k2 = 3*node[2], 
            k3 = 3*node[3],
            
            (* Compute elastic forces "F" on cell corners: *)
            restVol = cell.restVol,
            AI = cell.restInv,
            B = ShapeMatrix(pos, k0, k1, k2, k3),
            C = LR3x3.Mul(B, AI),
            C11 = C[0,0], C12 = C[0,1], C13 = C[0,2],
            C21 = C[1,0], C22 = C[1,1], C23 = C[1,2],
            C31 = C[2,0], C32 = C[2,1], C33 = C[2,2],
            (* T = (C^t)*C *)
            T11 = C11*C11 + C21*C21 + C31*C31,
            T12 = C11*C12 + C21*C22 + C31*C32,
            T13 = C11*C13 + C21*C23 + C31*C33,
            T22 = C12*C12 + C22*C22 + C32*C32,
            T23 = C12*C13 + C22*C23 + C32*C33,
            T33 = C13*C13 + C23*C23 + C33*C33,
            (* U = adj(T) *)
            U11 = T22*T33 - T23*T23,
            U12 = T12*T33 - T23*T13,
            U13 = T12*T23 - T22*T13,
            U22 = T11*T33 - T13*T13,
            U23 = T11*T23 - T12*T13,
            U33 = T11*T22 - T12*T12,
            Delta = T11*U11 - T12*U12 + T13*U13,
            Gamma = T11 + T22 + T33,
            dPdDelta =  restVol*cell.alpha/16.0D0*(Delta - 1.0D0/(Delta*Delta*Delta)),
            dPdSigma = -restVol*cell.beta/2.0D0,
            dPdGamma =  restVol*cell.beta*Gamma/3.0D0,
            (* d P/d T = (d P/d Delta)*(d Delta/d T) + 
                         (d P/d Sigma)*(d Sigma/d T) + 
                         (d P/d Gamma)*(d Gamma/d T) *)
            dPdT11 = 2.0D0*( dPdDelta*U11 + dPdSigma*(T22 + T33) + dPdGamma),
            dPdT22 = 2.0D0*( dPdDelta*U22 + dPdSigma*(T11 + T33) + dPdGamma),
            dPdT33 = 2.0D0*( dPdDelta*U33 + dPdSigma*(T11 + T22) + dPdGamma),
            dPdT12 = 2.0D0*(-dPdDelta*U12 - dPdSigma*T12),
            dPdT13 = 2.0D0*( dPdDelta*U13 - dPdSigma*T13),
            dPdT23 = 2.0D0*(-dPdDelta*U23 - dPdSigma*T23),
            (* obs.: the factor 2 in each dPdT comes from the dPdC_{ij}`s below *)
            (* d P/d C = (d P/d T)*(d T/d C) *)
            dPdC11 = C11*dPdT11 + C12*dPdT12 + C13*dPdT13,
            dPdC12 = C11*dPdT12 + C12*dPdT22 + C13*dPdT23,
            dPdC13 = C11*dPdT13 + C12*dPdT23 + C13*dPdT33,
            dPdC21 = C21*dPdT11 + C22*dPdT12 + C23*dPdT13,
            dPdC22 = C21*dPdT12 + C22*dPdT22 + C23*dPdT23,
            dPdC23 = C21*dPdT13 + C22*dPdT23 + C23*dPdT33,
            dPdC31 = C31*dPdT11 + C32*dPdT12 + C33*dPdT13,
            dPdC32 = C31*dPdT12 + C32*dPdT22 + C33*dPdT23,
            dPdC33 = C31*dPdT13 + C32*dPdT23 + C33*dPdT33,
            (* F = (d P/d C)*(d C/d B)*(d B/d q) *)
            Fx1 = dPdC11*AI[0,0] + dPdC12*AI[0,1] + dPdC13*AI[0,2],
            Fy1 = dPdC21*AI[0,0] + dPdC22*AI[0,1] + dPdC23*AI[0,2],
            Fz1 = dPdC31*AI[0,0] + dPdC32*AI[0,1] + dPdC33*AI[0,2],
            Fx2 = dPdC11*AI[1,0] + dPdC12*AI[1,1] + dPdC13*AI[1,2],
            Fy2 = dPdC21*AI[1,0] + dPdC22*AI[1,1] + dPdC23*AI[1,2],
            Fz2 = dPdC31*AI[1,0] + dPdC32*AI[1,1] + dPdC33*AI[1,2],
            Fx3 = dPdC11*AI[2,0] + dPdC12*AI[2,1] + dPdC13*AI[2,2],
            Fy3 = dPdC21*AI[2,0] + dPdC22*AI[2,1] + dPdC23*AI[2,2],
            Fz3 = dPdC31*AI[2,0] + dPdC32*AI[2,1] + dPdC33*AI[2,2],
            Fx0 = -(Fx1 + Fx2 + Fx3),
            Fy0 = -(Fy1 + Fy2 + Fy3),
            Fz0 = -(Fz1 + Fz2 + Fz3),

            (* Compute viscous forces "G" on cell corners: *)
            (* Suffix "D" denotes derivative w.r.t. time *)
            currVol = LR3x3.Det(B)/6.0D0,
            BD = ShapeMatrix(vel, k0, k1, k2, k3),
            CD = LR3x3.Mul(BD, AI),
            CD11 = CD[0,0], CD12 = CD[0,1], CD13 = CD[0,2],
            CD21 = CD[1,0], CD22 = CD[1,1], CD23 = CD[1,2],
            CD31 = CD[2,0], CD32 = CD[2,1], CD33 = CD[2,2],
            (* dT = transpose(CD)*C + transpose(C)*CD *)
            TD11 = 2.0D0*(CD11*C11 + CD21*C21 + CD31*C31),
            TD12 = CD11*C12 + CD21*C22 + CD31*C32 +
                   C11*CD12 + C21*CD22 + C31*CD32,
            TD13 = CD11*C13 + CD21*C23 + CD31*C33 +
                   C11*CD13 + C21*CD23 + C31*CD33,
            TD22 = 2.0D0*(CD12*C12 + CD22*C22 + CD32*C32),
            TD23 = CD12*C13 + CD22*C23 + CD32*C33 +
                   C12*CD13 + C22*CD23 + C32*CD33,
            TD33 = 2.0D0*(CD13*C13 + CD23*C23 + CD33*C33),
            Pi = TD11 + TD22 + TD33,
            dWdPi = currVol*(cell.eta1 + 4.0D0/3.0D0*cell.eta2)*Pi,
            dWdXi = -currVol*2.0D0*cell.eta2,
            (* d W/d T = (d W/d Pi)*(d Pi/d T) + (d W/d Xi)*(d Xi/d T) *)
            dWdTD11 = 2.0D0*(dWdPi + dWdXi*(TD22 + TD33)),
            dWdTD22 = 2.0D0*(dWdPi + dWdXi*(TD11 + TD33)),
            dWdTD33 = 2.0D0*(dWdPi + dWdXi*(TD11 + TD22)),
            dWdTD12 = 2.0D0*(-dWdXi*TD12),
            dWdTD13 = 2.0D0*(-dWdXi*TD13),
            dWdTD23 = 2.0D0*(-dWdXi*TD23),
            (* d W/d C = (d W/d T)*(d T/d C) *)
            dWdCD11 = C11*dWdTD11 + C12*dWdTD12 + C13*dWdTD13,
            dWdCD12 = C11*dWdTD12 + C12*dWdTD22 + C13*dWdTD23,
            dWdCD13 = C11*dWdTD13 + C12*dWdTD23 + C13*dWdTD33,
            dWdCD21 = C21*dWdTD11 + C22*dWdTD12 + C23*dWdTD13,
            dWdCD22 = C21*dWdTD12 + C22*dWdTD22 + C23*dWdTD23,
            dWdCD23 = C21*dWdTD13 + C22*dWdTD23 + C23*dWdTD33,
            dWdCD31 = C31*dWdTD11 + C32*dWdTD12 + C33*dWdTD13,
            dWdCD32 = C31*dWdTD12 + C32*dWdTD22 + C33*dWdTD23,
            dWdCD33 = C31*dWdTD13 + C32*dWdTD23 + C33*dWdTD33,
            (* G = (d W/d CD)*(d C/d B)*(d B/d q) *)
            Gx1 = dWdCD11*AI[0,0] + dWdCD12*AI[0,1] + dWdCD13*AI[0,2],
            Gy1 = dWdCD21*AI[0,0] + dWdCD22*AI[0,1] + dWdCD23*AI[0,2],
            Gz1 = dWdCD31*AI[0,0] + dWdCD32*AI[0,1] + dWdCD33*AI[0,2],
            Gx2 = dWdCD11*AI[1,0] + dWdCD12*AI[1,1] + dWdCD13*AI[1,2],
            Gy2 = dWdCD21*AI[1,0] + dWdCD22*AI[1,1] + dWdCD23*AI[1,2],
            Gz2 = dWdCD31*AI[1,0] + dWdCD32*AI[1,1] + dWdCD33*AI[1,2],
            Gx3 = dWdCD11*AI[2,0] + dWdCD12*AI[2,1] + dWdCD13*AI[2,2],
            Gy3 = dWdCD21*AI[2,0] + dWdCD22*AI[2,1] + dWdCD23*AI[2,2],
            Gz3 = dWdCD31*AI[2,0] + dWdCD32*AI[2,1] + dWdCD33*AI[2,2],
            Gx0 = -(Gx1 + Gx2 + Gx3),
            Gy0 = -(Gy1 + Gy2 + Gy3),
            Gz0 = -(Gz1 + Gz2 + Gz3)
          DO
            (* Add them into the foce vector: *)
            f[k0+0] := f[k0+0] - Fx0 - Gx0;
            f[k0+1] := f[k0+1] - Fy0 - Gy0;
            f[k0+2] := f[k0+2] - Fz0 - Gz0;
            f[k1+0] := f[k1+0] - Fx1 - Gx1;
            f[k1+1] := f[k1+1] - Fy1 - Gy1;
            f[k1+2] := f[k1+2] - Fz1 - Gz1;
            f[k2+0] := f[k2+0] - Fx2 - Gx2;
            f[k2+1] := f[k2+1] - Fy2 - Gy2;
            f[k2+2] := f[k2+2] - Fz2 - Gz2;
            f[k3+0] := f[k3+0] - Fx3 - Gx3;
            f[k3+1] := f[k3+1] - Fy3 - Gy3;
            f[k3+2] := f[k3+2] - Fz3 - Gz3;
          END
        END
      END
    END
  END ComputeInternalForces;

PROCEDURE ShapeMatrix(
    READONLY c: Coords;
    k0, k1, k2, k3: CARDINAL;
  ): LR3x3.T =
  (*
    
    If "q" is the system's position vector, computes the "B"
    matrix for some cell (see RLWL's thesis).  If "q" is the velocity
    vector, computes the time derivative of "B".
    
    The parameters "k0..k3" are coordinate indices (node indices × 3) 
    of the cell corners, in the canonical order. *)

  BEGIN
    WITH
      x0 = c[k0], y0 = c[k0+1], z0 = c[k0+2],
      x1 = c[k1], y1 = c[k1+1], z1 = c[k1+2],
      x2 = c[k2], y2 = c[k2+1], z2 = c[k2+2],
      x3 = c[k3], y3 = c[k3+1], z3 = c[k3+2]
    DO
      RETURN LR3x3.T{
        LR3.T{x1 - x0, y1 - y0, z1 - z0},
        LR3.T{x2 - x0, y2 - y0, z2 - z0},
        LR3.T{x3 - x0, y3 - y0, z3 - z0}
      }
    END
  END ShapeMatrix;
  
s.g := g;
      s.Kspring := Kspring;
      
      PROCEDURE TreatFixedNodes(s: T; v: Vector) =
  BEGIN
    IF s.nfixed > 0 THEN
      FOR i := 0 TO s.n-1 DO
        IF s.fixed[i] THEN
          WITH p = 3*i DO
            v[p] := 0.0D0; v[p+1] := 0.0D0; v[p+2] := 0.0D0
          END
        END
      END
    END
  END TreatFixedNodes;

BEGIN
END SystemDynamics.
