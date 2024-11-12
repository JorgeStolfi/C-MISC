INTERFACE Util;

(* Miscellaneous hacks *)

IMPORT Wr;

PROCEDURE Message(msg: TEXT);
  (*
    Writes "msg" to "stderr". *)

PROCEDURE Error(msg: TEXT);
  (*
    Writes "msg" to "stderr" and bombs out. *)

PROCEDURE Digits(n: CARDINAL): CARDINAL;
  (*
    Number of digits in "Fmt.Int(n)" *)

PROCEDURE ResetCPUTime();
PROCEDURE WriteCPUTime(wr: Wr.T);
  (*
    WriteCPUTime writes the CPU time elapsed since the 
    last call to "ResetCPUTime" or "WriteCPUTime", or since program
    startup if there were no such calls. *)

END Util.
