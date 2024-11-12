        FOR n := 0 TO s.nNodes-1 DO
          WITH fzn = f[3*i+2] DO
            fzn := fzn - s.nodeWeight[n]
          END
        END;
      gravity: LONGREAL;        (* Gravity. *)

      (* Node weights, for gravity: *)
      s.nodeWt := ComputeNodeWeights(s, gravity);
      
