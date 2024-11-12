(* 
Bipartite graph generator
Author: Jorge Stolfi
*)
MODULE JCSGraphs EXPORTS Main;

IMPORT Wr, Thread, Fmt, Math, Random;
IMPORT ParseParams, Process, Text;
FROM Stdio IMPORT stdout, stderr;

TYPE
  Kind = {Band, Hexa, Zipf, Fuzz, Grid, Worm, Rope, Puff};
  
  Options = RECORD
      NU, NV: CARDINAL;
      NE: CARDINAL;
      seed: CARDINAL;
      kind: Kind;
      dontScramble: BOOLEAN;
      reverse: BOOLEAN;
      comment: TEXT;
    END;
  
  Coins = Random.Default OBJECT METHODS
      init(coinSkip: CARDINAL): Coins := CoinsInit
    END;
    
  Edge = RECORD u, v: CARDINAL END;
  
PROCEDURE DoIt() =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      o = GetOptions(),
      e = NEW(REF ARRAY OF Edge, o.NE)^,
      coinSkip = o.seed * (o.NU + o.NV + o.NE)
    DO
(*      <* ASSERT o.NE <= o.NU * o.NV *> *)
      Wr.PutText(stderr, "Generating edges...\n");
      CASE o.kind OF 
      | Kind.Band => MakeBandGraph(o.NU, o.NV, coinSkip, e);
      | Kind.Hexa => MakeHexaGraph(o.NU, o.NV, coinSkip, e);
      | Kind.Zipf => MakeZipfGraph(o.NU, o.NV, coinSkip, e);
      | Kind.Fuzz => MakeFuzzGraph(o.NU, o.NV, coinSkip, e);
      | Kind.Grid => MakeGridGraph(o.NU, o.NV, coinSkip, e);
      | Kind.Worm => MakeWormGraph(o.NU, o.NV, coinSkip, e);
      | Kind.Rope => MakeRopeGraph(o.NU, o.NV, coinSkip, e);
      | Kind.Puff => MakePuffGraph(o.NU, o.NV, coinSkip, e);
      END;
      IF NOT o.dontScramble THEN
        Wr.PutText(stderr, "Relabeling vertices...\n");
        IF o.reverse THEN
          ReverseVertices(o.NU, o.NV, e)
        ELSE
          PermuteGraph(o.NU, o.NV, coinSkip, e)
        END
      END;    
      Wr.PutText(stderr, "Writing graph...\n");
      WriteGraph(stdout, o.comment, o.NU, o.NV, e);
      Wr.PutText(stderr, "Done.\n");
    END;
  END DoIt;
  
PROCEDURE MakeBandGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  (*
    This graph connects each vertex "u" in "U" to  
    the corresponding vertex "v" of "V", and possibly to other
    vertices of "V" in a band centered on "v",
    with probability decreasing with distance from "v".
  *)
  VAR nEdges: CARDINAL := 0;
      band: INTEGER;
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      coins = NEW(Coins).init(coinSkip := coinSkip),
      NE = NUMBER(e),
      (* Band width: *)
      dens = FLOAT(NE)/FLOAT(NU)/FLOAT(NV),
      BW = CEILING(FLOAT(NV)*dens*(2.0 - dens)),
      (* Coefficients for number of edges in band: *)
      fNE = FLOAT(NE, LONGREAL),
      fBW = FLOAT(BW, LONGREAL),
      A = FLOAT(NU, LONGREAL),
      B = 2.0d0*(A*fBW - fNE)/(fBW*(fBW - 1.0d0)),
      (* Denominator for cumulative edge density: *)
      FDen = fBW *(A - B*(fBW - 1.0d0)/2.0d0)
    DO
      IF NE = 0 THEN RETURN END;
      Wr.PutText(stderr, "BW = " & Fmt.Int(BW) & "\n");
      Wr.PutText(stderr, "A = " & Fmt.LongReal(A) & "\n");
      Wr.PutText(stderr, "B = " & Fmt.LongReal(B) & "\n");
       <* ASSERT BW > 0 *>
      <* ASSERT BW <= NV *>
      <* ASSERT NE <= NU * BW *>
      <* ASSERT B >= 0.0d0 *>
      FOR k := 0 TO BW-1 DO
        IF k MOD 2 = 0 THEN band := k DIV 2 ELSE band := - (k + 1) DIV 2 END;
        WITH
          (* Numerator for cumulative edge density: *)
          fk = FLOAT(k, LONGREAL),
          FNum = (fk + 1.0d0)*(A - B*fk/2.0d0),
          (* Ideal number of edges up to and including the "k"th band: *)
          iEdges = ROUND(FLOAT(NE, LONGREAL)*FNum/FDen),
          (* Number of edges to put in the "k"th band: *)
          d = MIN(NU, iEdges - nEdges),
          (* Jitter magnitude: *)
          Jitter = FLOAT(NU-d)/FLOAT(NU)
        DO
          FOR i := 0 TO d-1 DO 
            WITH
              u = ROUND(FLOAT(NU)*FLOAT(i)/FLOAT(d) + Jitter*RR(coins)) MOD NU,
              v = (ROUND(FLOAT(NV)*FLOAT(u)/FLOAT(NU)) + k) MOD NV
            DO
              PutEdge(u, v, e, nEdges)
            END
          END
        END
      END;
      <* ASSERT nEdges = NE *>
    END 
  END MakeBandGraph;
  
PROCEDURE MakeHexaGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  (*
    This graph consists of the union of edge-disjoint 
    random hexagons with shared vertices.
  *)
  VAR u, v: ARRAY [0..2] OF CARDINAL;
      nEdges: CARDINAL := 0;
  CONST PMax = 0.1;  (* Max prob of an hexagon sharing 2 vertices with some other *)
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      coins = NEW(Coins).init(coinSkip := coinSkip),
      NE = NUMBER(e),
      NH = NE DIV 6,
      blockSize = MIN(
        SQRT(FLOAT(NU)*FLOAT(NV)/FLOAT(NH)),
        MAX(3.0, 1.0 + (3.0/PMax)*FLOAT(NE)*(1.0/FLOAT(NU) + 1.0/FLOAT(NV)))
      ),
      NBU = CEILING(FLOAT(NU)/blockSize),
      NBV = CEILING(FLOAT(NV)/blockSize)
    DO
(*      <* ASSERT NE <= 6*NBU*NBV *>*)
      Wr.PutText(stderr, "blockSize = " & Fmt.Real(blockSize) & "\n");
      Wr.PutText(stderr, "NBU = " & Fmt.Int(NBU) & "\n");
      Wr.PutText(stderr, "NBV = " & Fmt.Int(NBV) & "\n");
      FOR ub := 0 TO NBU-1 DO 
        WITH
          (* Range of U block: *)
          uStart = FLOOR(FLOAT(NU)*FLOAT(ub)/FLOAT(NBU)),
          uStop = FLOOR(FLOAT(NU)*FLOAT(ub+1)/FLOAT(NBU)),
          (* Ideal number of edges up to and including this U block: *)
          iEdges = ROUND(FLOAT(NE)*FLOAT(ub+1)/FLOAT(NBU)),
          (* Number of edges to attach to this U block: *)
          de = iEdges - nEdges,
          (* Number of V blocks to connect to: *)
          db = (de + 5) DIV 6
        DO
          <* ASSERT db <= NBV *>
          <* ASSERT uStop - uStart >= 3 *>
          FOR kvb := 0 TO db - 1 DO
            WITH
              (* Which V block to connect to: *)
              vb = (FLOOR(FLOAT(NBV)*FLOAT(ub)/FLOAT(NBU)) + kvb) MOD NBV,
              (* Range of that V block: *)
              vStart = FLOOR(FLOAT(NV)*FLOAT(vb)/FLOAT(NBV)),
              vStop = FLOOR(FLOAT(NV)*FLOAT(vb+1)/FLOAT(NBV))
            DO
              <* ASSERT vStop - vStart >= 3 *>
              (* Generate a random bipartite hexagon between U and V blocks: *)
              RP(coins, 3, uStart, uStop - 1, u);
              RP(coins, 3, vStart, vStop - 1, v);
              <* ASSERT u[0] # u[1] AND u[0] # u[2] AND u[1] # u[2] *>
              <* ASSERT v[0] # v[1] AND v[0] # v[2] AND v[1] # v[2] *>
              FOR k := 0 TO 2 DO 
                IF nEdges < iEdges THEN PutEdge(u[k], v[(k+1)MOD 3], e, nEdges) END;
                IF nEdges < iEdges THEN PutEdge(u[k], v[(k+2)MOD 3], e, nEdges) END;
              END
            END
          END
        END 
      END;
      <* ASSERT nEdges = NE *>
    END 
  END MakeHexaGraph;
  
PROCEDURE MakeZipfGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  (*
    In this graph, the probability of an edge between
    vertex "u" and vertex "v" is roughly proportional 
    to 1/(u*v).  That is, most edges connect
    two low vertices, some connect a low vertex to a high vertex,
    and only a few connect two high vertices.
  *)
  VAR nEdges: CARDINAL := 0;
      vNext: CARDINAL;
  BEGIN
    WITH
      coins = NEW(Coins).init(coinSkip := coinSkip),
      LNU = LN(FLOAT(NU)),
      LNV = LN(FLOAT(NV+1)),
      NE = NUMBER(e)
    DO
      FOR u := 0 TO NU-1 DO
        WITH
          (* Ideal number of edges up to and including this U vertex: *)
          iEdges = ROUND(FLOAT(NE)*LN(FLOAT(u+1))/LNU),
          (* Degree of this U vertex: *)
          du = MIN(NV, iEdges - nEdges)
        DO
          vNext := 0;
          FOR k := 0 TO du-1 DO
            WITH
              vLo = MAX(vNext, ROUND(EX(LNV * FLOAT(k)/FLOAT(du))) - 1),
              vHi = MAX(vLo + 1, ROUND(EX(LNV * FLOAT(k+1)/FLOAT(du))) - 1),
              v = coins.integer(vLo, vHi-1)
            DO
              PutEdge(u, v, e, nEdges);
              vNext := v+1
            END
          END
        END
      END;
      <* ASSERT nEdges = NE *>
    END
  END MakeZipfGraph;
  
PROCEDURE LN(x: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.log(FLOAT(x, LONGREAL)))
  END LN;
  
PROCEDURE EX(x: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.exp(FLOAT(x, LONGREAL)))
  END EX;
  
PROCEDURE MakeFuzzGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  (*
    A "Fuzz" graph is a "Band" graph "G" with 
    a bunch leaves (unit-degree nodes) attached to 
    the vertices of G.  For low edge densities, 
    the number of leaves equals the number of non-leaves,
    so a perfect matching must use only the leaf edges.
  *)
  VAR nEdges: CARDINAL := 0;
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      NE = NUMBER(e),
      (* Maximum number of fuzz edges: *)
      B = FLOAT(2 - NU - NV, LONGREAL),
      C = FLOAT(NU, LONGREAL)*FLOAT(NV, LONGREAL) - FLOAT(NE, LONGREAL),
      MaxNF = MAX(0, FLOOR((-B + Math.sqrt(B*B - 4.0d0*C))/2.0d0 * 0.999999d0)),
      NF = MIN(MaxNF, MIN(NU DIV 2, NV DIV 2)),  (* Fuzz edges per side *)
      KU = NU - NF,     (* Non-leaf U vertices *)
      KV = NV - NF,     (* Non-leaf V vertices *)
      KE = NE - 2 * NF  (* Non-fuzz edges *)
    DO
      Wr.PutText(stderr, "KE = " & Fmt.Int(KE) & "\n");
      Wr.PutText(stderr, "KU = " & Fmt.Int(KU) & "\n");
      Wr.PutText(stderr, "KV = " & Fmt.Int(KV) & "\n");
      Wr.PutText(stderr, "NF = " & Fmt.Int(NF) & "\n");
      (* Ensure that there are not too many edges in kernel: *)
      IF  FLOAT(KE) > FLOAT(KU)*FLOAT(KV) THEN
        Wr.PutText(stderr, "Degree too high or too few vertices.\n");
        <* ASSERT FALSE *>
      END;
      (* Create non-fuzz edges: *)
      MakeBandGraph(KU, KV, coinSkip, SUBARRAY(e, 0, KE));
      nEdges := KE;
      (* Create fuzz edges: *)
      FOR k := 0 TO NF-1 DO
        WITH
          uf = KU + k,
          vf = KV + k,
          uo = FLOOR(FLOAT(KU)*FLOAT(k)/FLOAT(NF)),
          vo = FLOOR(FLOAT(KV)*FLOAT(k)/FLOAT(NF))
        DO
          PutEdge(uo, vf, e, nEdges);
          PutEdge(uf, vo, e, nEdges)
        END
      END;
      <* ASSERT nEdges = NE *>
    END 
  END MakeFuzzGraph;
  
PROCEDURE MakeGridGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  (*
    This graph is an approximate D-dimensional
    grid.  Each vertex "u" in "U" to 
    vertices "u+1", "u-1", "u+a", "u-a", "u+b", "u-b", ...
    where "1, a, b, ..." is a geometric progression.
    About 1/10 of these edges are omitted.
  *)
  VAR nEdges: CARDINAL := 0;
  CONST Jitter = 0.1011271243;
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      coins = NEW(Coins).init(coinSkip := coinSkip),
      NE = NUMBER(e),
      (* Max edges per dimension: *)
      MaxED = 2*FLOOR(FLOAT(MAX(NU,NV))/(1.0 + Jitter)*0.999999),
      (* Number of dimensions: *)
      ND = CEILING(FLOAT(NE)/FLOAT(MaxED) * 1.000001)
    DO
      Wr.PutText(stderr, "ND = " & Fmt.Int(ND) & "\n");
      (* Ensure that edges of different dimensions are distinct: *)
      <* ASSERT EX(LN(FLOAT(MIN(NU,NV)-1))/FLOAT(ND)) >= 2.0 *>
      FOR d := 0 TO ND-1 DO
        WITH
          (* Stride along dimension "d": *)
          uStep = EX(LN(FLOAT(NU-1))*FLOAT(d)/FLOAT(ND)) * 1.000001,
          vStep = EX(LN(FLOAT(NV-1))*FLOAT(d)/FLOAT(ND)) * 1.000001,
          (* Ideal number of edges up to and including dimension "d": *)
          iEdges = ROUND(FLOAT(NE)*FLOAT(d+1)/FLOAT(ND)),
          (* Number of edge pairs in this dimension: *)
          nep = (iEdges - nEdges + 1) DIV 2
        DO
          Wr.PutText(stderr, "  d = " & Fmt.Int(d) & "\n");
          Wr.PutText(stderr, "  uStep = " & Fmt.Real(uStep) & "\n");
          Wr.PutText(stderr, "  vStep = " & Fmt.Real(vStep) & "\n");
          Wr.PutText(stderr, "  nep = " & Fmt.Int(nep) & "\n");
          (* Ensure that the edges are not trivial: *)
          <* ASSERT uStep >= 1.0 *>
          <* ASSERT vStep >= 1.0 *>
          (* Ensure that the edges are all distinct: *)
          <* ASSERT FLOAT(nep)*(1.0 + Jitter) <= FLOAT(MAX(NU,NV)) *>
          FOR i := 0 TO nep-1 DO 
            WITH
              uf = FLOAT(NU)*FLOAT(i)/FLOAT(nep) + Jitter*RR(coins),
              vf = FLOAT(NV)*FLOAT(i)/FLOAT(nep) + Jitter*RR(coins),
              uo = ROUND(uf) MOD NU,
              vo = ROUND(vf) MOD NV,
              ud = ROUND(uf + uStep) MOD NU,
              vd = ROUND(vf + vStep) MOD NV
            DO
              IF nEdges < iEdges THEN PutEdge(uo, vd, e, nEdges) END;
              IF nEdges < iEdges THEN PutEdge(ud, vo, e, nEdges) END;
            END
          END
        END
      END;
      <* ASSERT nEdges = NE *>
    END 
  END MakeGridGraph;
  
PROCEDURE MakeWormGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  BEGIN
    WITH
      NE = NUMBER(e),
      (* Number of blocks on each side: *)
      K = ComputeWormLayers(NU, NV, NE)
    DO
      MakeGenWormGraph(NU, NV, K, coinSkip, e)
    END
  END MakeWormGraph;
  
PROCEDURE ComputeWormLayers(NU, NV, NE: CARDINAL): CARDINAL =
  VAR K: CARDINAL := 1;
      nLeaves: REAL := 1.0;
      M: CARDINAL := MIN(NU,NV)-1;
  BEGIN
    WITH
      avgDegree = FLOAT(2*NE)/FLOAT(NU+NV),
      branchRatio = avgDegree - 1.0
    DO
      IF branchRatio >= 2.0 THEN
        WHILE nLeaves < FLOAT(M) DO
          nLeaves := nLeaves * branchRatio;
          INC(K);
          M := (MIN(NU,NV)-1) DIV K
        END
      END
    END;
    WHILE K > 3 AND (3*K-2)*M > NE DO
      DEC(K);
      M := (MIN(NU,NV)-1) DIV K
    END;
    RETURN K
  END ComputeWormLayers;
  
PROCEDURE MakeRopeGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  BEGIN
    WITH
      NE = NUMBER(e),
      (* Number of blocks on each side: *)
      K = ComputeRopeLayers(NU, NV, NE)
    DO
      MakeGenWormGraph(NU, NV, K, coinSkip, e)
    END
  END MakeRopeGraph;

PROCEDURE ComputeRopeLayers(NU, NV, NE: CARDINAL): CARDINAL =
  VAR N: CARDINAL;
      K: CARDINAL;
      M: CARDINAL;
  BEGIN
    IF NU = NV THEN N := NU ELSE N := MIN(NU, NV)-1 END;
    K := N;
    M := N DIV K;
    WHILE K > 3 AND (K + (K-1)*M)*M < NE DO DEC(K); M := N DIV K END;
    RETURN K
  END ComputeRopeLayers;
  
PROCEDURE MakePuffGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  BEGIN
    WITH
      (* Number of blocks on each side: *)
      K = 3
    DO
      MakeGenWormGraph(NU, NV, K, coinSkip, e)
    END
  END MakePuffGraph;

PROCEDURE MakeGenWormGraph(
    NU, NV: CARDINAL;
    K: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  (*
    For low to moderate densities, most vertices of this
    graph are divided into 2*K layers of M vertices each,
    alternating between "U" and "V".  There are
    M edges between layers 0--1, 2--3, 4--5,.. (2K-2)--(2K-1)
    which must all be used in a maximal matching.
    There are also at most M^2 random "decoy" edges between layers 
    1--2, 3--4,.. (2K-3)--(2K-2).
  *)
  VAR nEdges: CARDINAL := 0;
      mEdges, dEdges: CARDINAL := 0; (* Matching and decoy edges *)

  PROCEDURE MakeMatching(uStart, vStart, n: CARDINAL) =
    (* 
      Creates a perfect parallel matching between 
      uStart + [0..n-1] and vStart + [0..n-1],
      updating "mEdges"
    *)
    BEGIN
      FOR i := 0 TO n-1 DO
        WITH u = uStart + i, v = vStart + i DO            
          PutEdge(u, v, e, nEdges);
          INC(mEdges)
        END
      END
    END MakeMatching;
    
  PROCEDURE MakeDecoys(coins: Coins; uStart, vStart, n, d: CARDINAL) =
    (* 
      Creates a jittered subgraph with "d" edges between vertices
      uStart + [0..n-1] and vStart + [0..n-1], updating dEdges.
    *)
    BEGIN
      IF d = 0 THEN RETURN END;
      WITH
        Jitter = FLOAT(n*n-d)/FLOAT(n*n)
      DO
        FOR k := 0 TO d-1 DO
          WITH
            uv = FLOOR(FLOAT(n*n)*FLOAT(k)/FLOAT(d) + Jitter*RR(coins)),
            u = uStart + uv DIV n,
            v = vStart + uv MOD n
          DO
            PutEdge(u, v, e, nEdges);
            INC(dEdges)
          END
        END
      END
    END MakeDecoys;

  PROCEDURE MakeOverflow(coins: Coins; uStart, vStart, nU, nV, d: CARDINAL) =
    (* 
      Creates a jittered subgraph with "d" edges between vertices
      uStart + [0..nU-1] and vStart + [0..nV-1].
    *)
    BEGIN
      IF d = 0 THEN RETURN END;
      <* ASSERT nU*nV > 0 *>
      WITH
        Jitter = FLOAT(nU*nV-d)/FLOAT(nU*nV)
      DO
        FOR k := 0 TO d-1 DO
          WITH
            uv = FLOOR(FLOAT(nU*nV)*FLOAT(k)/FLOAT(d) + Jitter*RR(coins)),
            u = uStart + uv DIV nV,
            v = vStart + uv MOD nV
          DO
            PutEdge(u, v, e, nEdges)
          END
        END
      END
    END MakeOverflow;

  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      coins = NEW(Coins).init(coinSkip := coinSkip),
      NE = NUMBER(e),
      (* Size of each vertex block: *)
      M = ComputeWormThickness(K, NU, NV, NE),
      (* Number of decoy edges in worm itself: *)
      ND = CEILING(MIN(FLOAT(K-1)*FLOAT(M)*FLOAT(M), FLOAT(NE) - FLOAT(K)*FLOAT(M)))
    DO
      <* ASSERT M >= 0 *>
      <* ASSERT ND >= 0 *>
      Wr.PutText(stderr, "K = " & Fmt.Int(K) & "\n");
      Wr.PutText(stderr, "M = " & Fmt.Int(M) & "\n");
      Wr.PutText(stderr, "ND = " & Fmt.Int(ND) & "\n");
      
      (* Build worm: *)
      FOR b := 0 TO K-1 DO
        WITH
          (* NS = 0 for middle layers, = 1 for other layers *)
          NS = 1 - (b + 1) DIV K
        DO
          FOR side := 0 TO NS DO
            IF b MOD 2 = 0 THEN
              (* Even blocks: matching with next block *)
              MakeMatching((b+side)*M, (b+NS-side)*M, M)
            ELSE
              (* Odd blocks: decoy edges *)
              WITH
                t = b + side,
                DE = ROUND(FLOAT(ND)*FLOAT(t)/FLOAT(K-1)) - dEdges
              DO
                <* ASSERT FLOAT(DE) <= FLOAT(M)*FLOAT(M) *>
                MakeDecoys(coins, (b+side)*M, (b+NS-side)*M, M, DE)
              END
            END
          END
        END
      END;
      <* ASSERT mEdges = K*M *>
      
      (* Add overflow edges: *)
      MakeOverflow(coins, K*M, K*M, NU-K*M, NV-K*M, NE - nEdges);
      Wr.PutText(stderr,"OverFlow Edges = " & Fmt.Int(NE - nEdges) & "\n");
      <* ASSERT nEdges = NE *>
    END 
  END MakeGenWormGraph;
  
PROCEDURE ComputeWormThickness(K, NU, NV, NE: CARDINAL): CARDINAL =
  (*
    Computes the block size for a "Worm" graph with 2*K layers.
  *)
  VAR M: CARDINAL;
  BEGIN
    IF NU = NV THEN 
      M := NU DIV K
    ELSE
      (* Must leave at least one "overflow" vertex on each side *)
      M := (MIN(NU,NV)-1) DIV K
    END;
    IF M > 0 THEN
      WITH
        fK = FLOAT(K, LONGREAL),
        fM = FLOAT(M, LONGREAL),
        A = fK*(fK + 1.0d0) - 1.0d0,
        B = fK * FLOAT(1 - NU - NV, LONGREAL),
        C = FLOAT(NU, LONGREAL)*FLOAT(NV, LONGREAL) - FLOAT(NE, LONGREAL)
      DO
        IF (A*fM + B)*fM + C < 0.000001d0 THEN
          (* Must increase overflow area to accomodate excess edges: *)
          M := FLOOR((-B - Math.sqrt(B*B - 4.0d0 * A * C))/(2.0d0 * A)*0.999999d0);
        END
      END;
    END;
    M := MAX(0, M);
(*    <* ASSERT NE <= K*M + (K-1)*M*M + (NU-K*M)*(NV-K*M) *> *)
    RETURN M
  END ComputeWormThickness;
  
PROCEDURE PermuteGraph(
    NU, NV: CARDINAL;
    coinSkip: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  BEGIN
    WITH
      coins = NEW(Coins).init(coinSkip := coinSkip),
      NE = NUMBER(e),
      uLab = NEW(REF ARRAY OF CARDINAL, NU)^,
      vLab = NEW(REF ARRAY OF CARDINAL, NV)^
    DO
      FOR u := 0 TO NU-1 DO uLab[u] := u END; Scramble(coins, uLab);
      FOR v := 0 TO NV-1 DO vLab[v] := v END; Scramble(coins, vLab);
      (* Translate edges: *)
      FOR i := 0 TO NE-1 DO 
        WITH ei = e[i] DO
          ei.u := uLab[ei.u]; ei.v := vLab[ei.v]
        END;
      END;
    END;
  END PermuteGraph;
  
PROCEDURE ReverseVertices(
    NU, NV: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  BEGIN
    WITH
      NE = NUMBER(e)
    DO
      (* Translate edges: *)
      FOR i := 0 TO NE-1 DO 
        WITH ei = e[i] DO
          ei.u := NU - 1 - ei.u; 
          ei.v := NV - 1 - ei.v
        END;
      END;
    END;
  END ReverseVertices;
  
<*UNUSED*>
PROCEDURE PermuteDegrees(
    VAR d: ARRAY OF CARDINAL; 
    READONLY lab: ARRAY OF CARDINAL;
  ) =
  BEGIN
    WITH
      N = NUMBER(d),
      t = NEW(REF ARRAY OF CARDINAL, N)^
    DO
      <* ASSERT NUMBER(lab) = N *>
      FOR i := 0 TO N-1 DO t[lab[i]] := d[i] END;
      d := t
    END
  END PermuteDegrees;
  
PROCEDURE WriteGraph(
    wr: Wr.T;
    comment: TEXT;
    NU, NV: CARDINAL;
    VAR e: ARRAY OF Edge;
  ) =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    WITH
      NE = NUMBER(e),
      du = CEILING(Math.log10(FLOAT(NU, LONGREAL))),
      dv = CEILING(Math.log10(FLOAT(NV, LONGREAL)))
    DO
      Wr.PutText(wr, "p bipartite-graph");
      Wr.PutText(wr, " " & Fmt.Int(NU));
      Wr.PutText(wr, " " & Fmt.Int(NV));
      Wr.PutText(wr, " " & Fmt.Int(NE));
      Wr.PutText(wr, "\n");
      Wr.PutText(wr, "c");
      Wr.PutText(wr, " " & comment);
      Wr.PutText(wr, "\n");
      FOR i := 0 TO NE-1 DO 
        Wr.PutText(wr, "a");
        Wr.PutText(wr, " " & Fmt.Pad(Fmt.Int(e[i].u + 1), du));
        Wr.PutText(wr, " " & Fmt.Pad(Fmt.Int(e[i].v + 1), dv)); 
        Wr.PutText(wr, "\n");
        Wr.Flush(wr);
      END;
    END
  END WriteGraph;
  
<*UNUSED*>
PROCEDURE Zipf(N: CARDINAL): REF ARRAY OF REAL =
  (* Returns a sorted Zipf distribution *)
  VAR sum: REAL := 0.0;
  BEGIN
    WITH
      r = NEW(REF ARRAY OF REAL, N),
      p = r^
    DO
      FOR i := 0 TO N-1 DO p[i] := 1.0/FLOAT(i+1); sum := sum + p[i] END;
      FOR i := 0 TO N-1 DO p[i] := p[i]/sum END;
      RETURN r
    END
  END Zipf;

<*UNUSED*>
PROCEDURE RB(coins: Coins): BOOLEAN =
  BEGIN
    RETURN coins.boolean()
  END RB;

PROCEDURE RP(
    coins: Coins; 
    n: CARDINAL; 
    lo, hi: CARDINAL; 
    VAR v: ARRAY OF CARDINAL;
  ) =
  (*
    Puts in "v[0..n-1]" a random sequence of "n" distinct
    elements from [lo..hi].
  *)
  BEGIN
    RS(coins, n, lo, hi, v);
    Scramble(coins, v);
  END RP;
  
PROCEDURE RS(
    coins: Coins; 
    n: CARDINAL; 
    lo, hi: CARDINAL; 
    VAR v: ARRAY OF CARDINAL;
  ) =
  (*
    Puts in "v[0..n-1]" a random subset of "n" distinct
    elements from [lo..hi], in increasing order.
  *)
  VAR t, e: INTEGER;
  BEGIN
    <* ASSERT hi >= lo + n - 1 *>
    FOR k := 0 TO n-1 DO
      e := coins.integer(lo + k, hi);
      t := k;
      WHILE t > 0 AND v[t-1] >= e DO v[t] := v[t-1]; t := t-1; e := e - 1 END;
      v[t] := e
    END
  END RS;
  
PROCEDURE Scramble(coins: Coins; VAR v: ARRAY OF CARDINAL) =
  (* Permutes "v" in random order *)
  VAR t: CARDINAL;
  BEGIN
    FOR k := LAST(v) TO 1 BY -1 DO 
      WITH i = coins.integer(0, k) DO 
        IF i # k THEN
          t := v[i]; v[i] := v[k]; v[k] := t
        END
      END
    END
  END Scramble;
  
PROCEDURE SQRT(x: REAL): REAL =
  BEGIN
    RETURN FLOAT(Math.sqrt(FLOAT(x, LONGREAL)))
  END SQRT;

PROCEDURE PutEdge(
    u, v: CARDINAL;
    VAR e: ARRAY OF Edge;
    VAR nEdges: CARDINAL;
  ) =
  BEGIN
    e[nEdges] := Edge{u, v};
    INC(nEdges)
  END PutEdge;
  
PROCEDURE RR(coins: Coins): REAL = 
  BEGIN
    RETURN coins.real()
  END RR;

PROCEDURE GetOptions (): Options =
  <* FATAL Thread.Alerted, Wr.Failure *>
  VAR o: Options;
      kt: TEXT;
  CONST
    MaxVertices = 512*512;
    MaxEdges = 16*MaxVertices;
    MaxSeed = 999;
  BEGIN
    WITH pp = NEW(ParseParams.T).init(stderr) DO
      TRY
        pp.getKeyword("-kind");
        kt := pp.getNext();
        IF Text.Equal(kt, "Band") THEN
          o.kind := Kind.Band
        ELSIF Text.Equal(kt, "Fuzz") THEN
          o.kind := Kind.Fuzz
        ELSIF Text.Equal(kt, "Hexa") THEN
          o.kind := Kind.Hexa
        ELSIF Text.Equal(kt, "Zipf") THEN
          o.kind := Kind.Zipf
        ELSIF Text.Equal(kt, "Grid") THEN
          o.kind := Kind.Grid
        ELSIF Text.Equal(kt, "Worm") THEN
          o.kind := Kind.Worm
        ELSIF Text.Equal(kt, "Rope") THEN
          o.kind := Kind.Rope
        ELSIF Text.Equal(kt, "Puff") THEN
          o.kind := Kind.Puff
        ELSE
          pp.error("bad kind \"" & kt & "\"");
        END;

        o.comment := "kind: " & kt;

        pp.getKeyword("-NU");
        o.NU := pp.getNextInt(1, MaxVertices);

        pp.getKeyword("-NV");
        o.NV := pp.getNextInt(1, MaxVertices);

        pp.getKeyword("-NE");
        o.NE := pp.getNextInt(1, MaxEdges);

        pp.getKeyword("-seed");
        o.seed := pp.getNextInt(0, MaxSeed);

        o.comment := o.comment & " seed: " & Fmt.Int(o.seed);

        o.dontScramble := pp.keywordPresent("-dontScramble");
        IF NOT o.dontScramble THEN
          o.reverse := pp.keywordPresent("-reverse")
        ELSE
          o.reverse := FALSE
        END;

        pp.finish();                                       
      EXCEPT                                                            
      | ParseParams.Error =>                                              
          Wr.PutText(stderr, "Usage: JCSGraphs \\\n");
          Wr.PutText(stderr, "  -NU <num>  -NV <num>  -NE <num>  -seed <num>\\\n");
          Wr.PutText(stderr, "  [ -dontScramble | -reverse ]\\\n");
          Wr.PutText(stderr, "  -kind { Band | Hexa | Zipf \\\n");
          Wr.PutText(stderr, "        | Fuzz | Grid | Worm \\\n");
          Wr.PutText(stderr, "        | Rope | Puff }\n");
          Process.Exit (1);
      END;
    END;
    RETURN o
  END GetOptions;

PROCEDURE CoinsInit(coins: Coins; coinSkip: CARDINAL): Coins =
  (*
    This routine initializes the random number generator
    so that each seed generates a different but fixed sequence.  
    Unfortunately, the random-number generator in release 3.5
    does not provide a method to explicitly set the seed.
    So we merely discarding the first "coinSkip"
    elements of the Random.Default's fixed sequence. *)
  BEGIN
    WITH c = NARROW(coins, Random.Default).init(fixed := TRUE) DO
      FOR i := 1 TO coinSkip DO EVAL c.integer() END;
      RETURN c
    END
  END CoinsInit;

BEGIN
  DoIt()
END JCSGraphs.