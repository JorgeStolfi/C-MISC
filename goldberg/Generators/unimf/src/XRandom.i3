INTERFACE XRandom;

(* A variant of "Random.T" that can be initialized with a given seed. *)

(* Created June 1997 by J.Stolfi. *)
(* Based on "Random.i3" by Bill Kalsow. *)
(* See the copyright notice at the end of this file. *)

IMPORT Random;

TYPE
  T <: Public;
  Public = Random.T OBJECT METHODS
      init(seed: Seed): T;
    END;
  
  Seed = [0 .. (LAST(INTEGER) DIV 10) * 8 - 1];

END XRandom.

(* Copyright (C) 1997, Institute of Computing, UNICAMP. *)
