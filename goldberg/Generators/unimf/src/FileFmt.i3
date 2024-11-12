INTERFACE FileFmt;

(* Writing and parsing version-checked file headers, footers, and comments. *)
(* Created April 1997 by J. Stolfi  *)
(* See the copyright notice at the end of this file. *)

IMPORT Wr, Rd;

(* FILE HEADERS AND FOOTERS *)

(*
  These procedures help maintaining files with typed and versioned
  headers and footers, like this:
  
     begin mytable (format of 96-11-29)
     [...stuff...]
     end mytable
     
  The strings "mytable" and "96-11-29" are client-specified,
  and checked on input.  These begin/end pairs may be nested, etc.
  
  The procedures abort on any errors.

  The "ReadHeader" and "ReadFooter" procedures ignore whitespace
  before or after words, but require the entire header or footer to be
  contained in one line, and terminated by a newline.

*)

PROCEDURE WriteHeader(wr: Wr.T; type, version: TEXT);
PROCEDURE WriteFooter(wr: Wr.T; type: TEXT);

PROCEDURE ReadHeader(rd: Rd.T; type, version: TEXT);
PROCEDURE ReadFooter(rd: Rd.T; type: TEXT);

PROCEDURE MakeHeader(type, version: TEXT): TEXT;
PROCEDURE MakeFooter(type: TEXT): TEXT;

(* COMMENTS: *)

(*
  
  These routines write and parse a comment text, consisting of zero or
  more lines marked by a given "prefix" character:
  
    | Blah blah blah
    |   blah blah
    | and more blah.
  
  The "prefix" must be the first character of the line, and is
  normally followed by a blank. *)

PROCEDURE WriteComment(wr: Wr.T; comment: TEXT; prefix: CHAR);
  (* 
    Writes the given "comment" text to "wr", with a "prefix" character
    and a blank in front of every line.  Supplies a final '\n' if 
    the text is non-empty but does not end with newline. *)

PROCEDURE ReadComment(rd: Rd.T; prefix: CHAR): TEXT;
  (*
    Reads zero or more lines from "rd" that begin with the  "prefix"
    character, strips the leading "prefix" and the following blank
    (if present) from each line, and returns them as a single TEXT,
    with each line terminated by a newline. *)

END FileFmt.
