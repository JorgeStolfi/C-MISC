MODULE FileFmt;

(* Created April 1997 by J. Stolfi  *)
(* See the copyright notice at the end of this file. *)

IMPORT Rd, Wr, TextRd, TextWr, Lex;
FROM Thread IMPORT Alerted;
FROM Rd IMPORT EndOfFile;

CONST 
  Spaces = SET OF CHAR{' ', '\t'};

PROCEDURE WriteHeader(wr: Wr.T; type, version: TEXT) =
  <* FATAL Wr.Failure, Alerted *>
  BEGIN
    Wr.PutText(wr, "begin ");
    Wr.PutText(wr, type);
    Wr.PutText(wr, " (format of ");
    Wr.PutText(wr, version);
    Wr.PutText(wr, ")\n");
  END WriteHeader;

PROCEDURE WriteFooter(wr: Wr.T; type: TEXT) =
  <* FATAL Wr.Failure, Alerted *>
  BEGIN 
    Wr.PutText(wr, "end ");
    Wr.PutText(wr, type);
    Wr.PutText(wr, "\n") 
  END WriteFooter;

PROCEDURE MakeHeader(type, version: TEXT): TEXT = 
  BEGIN
    WITH wr = NEW(TextWr.T).init() DO
      WriteHeader(wr, type, version);
      RETURN TextWr.ToText(wr)
    END
  END MakeHeader;

PROCEDURE MakeFooter(type: TEXT): TEXT =
  BEGIN
    WITH wr = NEW(TextWr.T).init() DO
      WriteFooter(wr, type);
      RETURN TextWr.ToText(wr)
    END
  END MakeFooter;

PROCEDURE ReadHeader(rd: Rd.T; type, version: TEXT) =
  <* FATAL Rd.Failure, Alerted, Lex.Error *>
  CONST Spaces = SET OF CHAR{' ', '\t'};
  BEGIN
    Lex.Skip(rd, Lex.Blanks);
    Lex.Match(rd, "begin");
    Lex.Skip(rd, Spaces);
    Lex.Match(rd, type);
    Lex.Skip(rd, Spaces);
    Lex.Match(rd, "(format of");
    Lex.Skip(rd, Spaces);
    Lex.Match(rd, version);
    Lex.Skip(rd, Spaces);
    Lex.Match(rd, ")");
    Lex.Skip(rd, Spaces);
    Lex.Match(rd, "\n");
  END ReadHeader;

PROCEDURE ReadFooter(rd: Rd.T; type: TEXT) =
  <* FATAL Rd.Failure, Alerted, Lex.Error *>
  BEGIN 
    Lex.Skip(rd, Lex.Blanks);
    Lex.Match(rd, "end ");
    Lex.Skip(rd, Spaces);
    Lex.Match(rd, type);
    Lex.Skip(rd, Spaces);
    Lex.Match(rd, "\n");
  END ReadFooter;

EXCEPTION MissingFinalNewLine;

PROCEDURE WriteComment(wr: Wr.T; comment: TEXT; prefix: CHAR) =
  
  VAR rd: Rd.T := TextRd.New(comment);
  
  PROCEDURE CopyLine() RAISES {EndOfFile} =
  (*
    Copy one line from "rd" to "wr", prefixed by "prefix"
    and a space. Supplies a final '\n' if the next line exists but
    does not end with newline. Raises "EndOfFile" if there are no more
    lines in "rd". *)
    
    <* FATAL Rd.Failure, Wr.Failure, Alerted *>
    VAR c: CHAR;
    BEGIN
      c := Rd.GetChar(rd); (* If EOF here, propagate to caller *)
      Wr.PutChar(wr, prefix);
      Wr.PutChar(wr, ' ');
      Wr.PutChar(wr, c);
      WHILE c # '\n' DO
        TRY c := Rd.GetChar(rd) EXCEPT EndOfFile => c := '\n' END;
        Wr.PutChar(wr, c)
      END
    END CopyLine;

  BEGIN
    TRY LOOP CopyLine() END EXCEPT EndOfFile => (* Ok *) END;
  END WriteComment;

PROCEDURE ReadComment(rd: Rd.T; prefix: CHAR): TEXT =

  VAR wr: Wr.T := NEW(TextWr.T).init();

  PROCEDURE CopyLine() RAISES {EndOfFile} =
  (*
    Copy one comment line from "rd" to "wr", removing the "prefix"
    and the following blank, but leaving the final (mandatory) newline.
    Raises EndOfFile if "rd" is exhausted or the next char is 
    not "prefix". *)
    <* FATAL Rd.Failure, Wr.Failure, Alerted *>
    <* FATAL MissingFinalNewLine *>
    VAR c: CHAR;
    BEGIN
      c := Rd.GetChar(rd); (* If EOF here, propagate to caller *)
      IF c # prefix THEN Rd.UnGetChar(rd); RAISE EndOfFile END;
      TRY c := Rd.GetChar(rd) EXCEPT EndOfFile => RAISE MissingFinalNewLine END;
      IF c = ' ' THEN 
        TRY c := Rd.GetChar(rd) EXCEPT EndOfFile => RAISE MissingFinalNewLine END
      END;
      WHILE c # '\n' DO
        Wr.PutChar(wr, c);
        TRY c := Rd.GetChar(rd) EXCEPT EndOfFile => RAISE MissingFinalNewLine END;
      END;
      Wr.PutChar(wr, c);
    END CopyLine;

  BEGIN
    TRY LOOP CopyLine() END EXCEPT EndOfFile => (* Ok *) END;
    RETURN TextWr.ToText(wr)
  END ReadComment;

BEGIN
END FileFmt.
