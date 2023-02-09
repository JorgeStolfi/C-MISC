MODULE MaintainAutomaton EXPORTS Main;

(* Maintenance operations on wordlist automata. *)
(* See the copyright and disclaimer note at the end of this file. *)

(*
  This program performs maintenance operations on wordlist automata.

  Generally speaking, the program reads or builds an initial
  automaton, performs some modifications on it, prints some
  reports and auxiliary listings, and optionally writes
  the resulting automaton to a "dump" file.
  
  The initial automaton can be specified by giving the command:
  
    -load FILENAME
    
        Read the initial automaton from FILENAME, in the ".dmp" 
        format (see "Reduced.Dump").
  
  If no "-load" command is given, the initial automaton is empty.
  In any case, this initial automaton can be modified by 
  at most one of the following commands:

    -add FILENAME
    
        Read words from FILENAME, add them to the initial automaton.

    -remove FILENAME
    
        Read words from FILENAME, remove them from the initial automaton.

    -addRemove FILENAME
        
        Read words from FILENAME, add/remove them to/from the initial
        automaton.  Each word in FILENAME must be preceded by 
        '+' (add) or '-' (remove), and optional whitespace.

  Any of these operations may be combined with any subset of the following:

    -dump FILENAME
    
        Write the automaton to FILENAME, in ".dmp" format
        (see "Reduced.Dump").
     
    -check FILENAME
    
        Read words from FILENAME, and print on "stdout" those that are
        not accepted by the automaton.
               
    -spell FILENAME

        Write to FILENAME all words accepted by the automaton.

    -print FILENAME
        
        Write to FILENAME a readable description of the automaton.

    -pickle FILENAME
    
        Write to FILENAME a Modula-3 pickle of the automaton.

  These commands are applied to the modified automaton that resulted
  from the "-add"/"-remove"/"-addRemove" operations, or to the initial
  automaton if no modification was requested.

  The files read by "-add", "-remove", "-adRemove", and "-check" must
  contain one word per line.  Valid word characters are all printing
  ISO-Latin-1 characters ('\041' to '\176' and '\241' to '\377') plus
  the plain ASCII space ('\040'). Leading and trailing blanks are
  discarded, but embedded spaces are significant and treated as
  letters.  In particular, a blank line is treated as the empty word.

  The "-spell" command writes the words in the same format, without
  leading or trailing blanks, in lexicographic order.
  
  Other commands and options are:
  
    -help
    
        Print a summary of the command line arguments.
        (All other commnds are ignored.)
   
    -log FILENAME
    
        Redirect all error and informative messages (except command line 
        parsing errors) to FILENAME.  Default is "stderr".

    -stats          
        
        Print the final automaton's statistics to the "log" file.

    -comment TEXT     
    
        Include TEXT in the "-dump" file, as a comment.

    -commentFile FILENAME 
    
        Copy FILENAME into the "-dump" file, as comment.
    
    -dagSize <num>      
    
        Pre-allocate this many transitions for the initial automaton
        (default: 100 000 transitions).  The allocated area is
        automatically expanded as needed, anyway; so this optimization
        does not affect the result.
                         
    -reportStep NUM   
    
        For the "-add", "-remove", "-addRemove", and "-check"
        commands: print a status line after processing every NUM input
        words (default: every 10 000 words)
                           
    -noRedundant        
    
        For the "-add", "-remove", and "-addRemove" commands: do not
        issue warnings when processing "redundant" words (words to
        include which are already in the automaton, or words to remove
        which aren't there.)

HISTORY:

  13.apr.92: [TK] Modified for compatibility with version 2.04 of Modula-3 and the
             new interfaces.

  02.nov.92: [TK] "-noredundant" option added.   

  06.aug.95: [JS] converted to Modula-3 release 3.5.3.  Subtantial rewrite
             of command line processing.  Fixed comments to better match the
             code, etc.
             
  24.sep.95: [JS] Removed the "-build" command (use "-add" without "-load").
             Changed command syntax to require FILENAME for the commands
             "-check", "-add", "-remove", "-addRemove" (which formerly
             read from stdin) and "-spell" (which formerly wrote to stdout).
             
  01.nov.96  [JS] Added a quick fix to the word length histogram code
             so that words longer than "MaxWordLength" are counted too.
             
  19.jan-97  [JS] Removed the historical compatibility hacks ("-dumpin",
             "-dumpout", "-mess", etc.  Upgraded to "libm3reduced-2".

*)

IMPORT ParseParams, Process, Wr, Rd, Text, TextWr, Fmt;
IMPORT Pickle, Thread;
IMPORT Reduced, ParamUtil, Util, Encoding, PlainEncoding;
FROM Stdio IMPORT stdout, stderr;
FROM Text IMPORT Empty;
FROM Basics IMPORT BOOL, NAT, Abort, Done, String, Letter, WrNat;

CONST HelpText =
  "options:\n" &
  "  MaintainAutomaton \\\n" &
  "    [ -log FILE.log ] \\\n" &
  "    [ -load FILE.dmp ] \\\n" &
  "    [ -dagSize NNN ] \\\n" &
  "    [ -add FILE.dic | -remove FILE.dic | -addRemove FILE.fix  ] \\\n" &
  "    [ -check FILE.dic ] \\\n" &
  "    [ -noRedundant ] \\\n" &
  "    [ -reportStep NNN ] \\\n" &
  "    [ -stats ] \\\n" &
  "    [ -dump FILE.dmp [ -comment \"blah blah\" | -commentFile FILE.cmt ] ] \\\n" &
  "    [ -spell FILE.out ] \\\n" &
  "    [ -print FILE.prn ] \\\n" &
  "    [ -pickle FILE.pkl ]\n";

<* FATAL Thread.Alerted *>

CONST
  Version = "3.5.3-1";
  MaxWordLength = 100;
  DefaultDagSize = 100000;
  DefaultReportStep = 10000;
       
VAR
  encoding: Encoding.T := PlainEncoding.New();
  
  arg: RECORD
  
      log: Wr.T;         (* Log file, or "stderr" *)
      cmd: TEXT;         (* Formatted command line arguments *)

      load: TEXT;        (* Input ".dmp" filename, or "" if none *)
      dagSize: NAT;      (* Initial allocation size *)
      modify: TEXT;      (* Input filename for modify op, or "" if none *)
      op: ModifyOp;      (* Operation to apply, if "modify" is not "". *)
      noRedundant: BOOL; (* TRUE to ignore redundant add/remove lines *)
      reportStep: NAT;   (* Progress report period *)
      stats: BOOL;       (* TRUE to print automaton statistics to "arg.log" *)

      (* Output options: *)
      dump: TEXT;        (* Output ".dmp" filename, or "" if none *) 
      doc: TEXT;         (* New comment text *)
      check: TEXT;       (* Input filename, or "" if none *) 
      spell: TEXT;       (* Filename for final automaton wordlist, or "" if none *)
      print: TEXT;       (* Filename for final automaton printout, or "" if none *)
      pickle: TEXT;      (* Filename for final automaton pickle, or "" if none *)
    END;

TYPE
  ModifyOp = {
    Add,       (* Add words *)            
    Remove,    (* Remove words *)        
    AddRemove  (* Add/remove words, depending on first char *) 
  };

CONST
  ModifyOpFlag = ARRAY ModifyOp OF TEXT {
      (* ModifyOp.Add:       *) "-add", 
      (* ModifyOp.Remove:    *) "-remove", 
      (* ModifyOp.AddRemove: *) "-addRemove"
    };

  CONST ModifyOpMessage = ARRAY ModifyOp OF TEXT{
      (* Op.Add:       *) "adding words to the automaton",
      (* Op.Remove:    *) "removing words from the automaton",
      (* Op.AddRemove: *) "adding/removing words from the automaton"
    };

PROCEDURE DoIt() =
  VAR
    aut: Reduced.T;
    oldDoc: TEXT;
  BEGIN
    GetCommandLineArguments();
    
    StartedMessage("MaintainAutomaton " & Version & " " & Util.FmtDate(Util.GetDate()));
    
    DoShowArgs(arg.cmd);

    (* Create the automaton: *)
    IF NOT Empty(arg.load) THEN
      aut := DoLoad(Util.OpenRd(arg.load, "-load"), arg.dagSize);
      oldDoc := aut.doc;
    ELSE 
      oldDoc := "";
      aut := Reduced.New(arg.dagSize)
    END;
    
    IF NOT Empty(arg.modify) THEN
      WITH rd = Util.OpenRd(arg.modify, ModifyOpFlag[arg.op]) DO
        DoModify(aut, rd, arg.op, arg.noRedundant, arg.reportStep)
      END
    END;

    (* Write automaton statistics: *)
    IF arg.stats THEN 
      DoWriteStats(aut)
    END;

    IF NOT Empty(arg.check) THEN
      DoCheck(aut, Util.OpenRd(arg.check, "-check"), stdout, arg.reportStep);
    END;

    (* Write the ".dmp" automaton: *)
    IF NOT Empty(arg.dump) THEN
      aut.doc := MakeDoc(oldDoc, arg.cmd, arg.doc);
      DoDump(aut, Util.OpenWr(arg.dump, "-dump"));
    END;

    (* Write the word list: *)
    IF NOT Empty(arg.spell) THEN
      DoSpell(aut, Util.OpenWr(arg.spell, "-spell"), arg.reportStep)
    END;

    (* Write the readable automaton: *)
    IF NOT Empty(arg.print) THEN
      DoPrint(aut, Util.OpenWr(arg.print, "-print"));
    END;
    
    (* Write the pickled automaton: *)
    IF NOT Empty(arg.pickle) THEN
      DoPickle(aut, Util.OpenWr(arg.pickle, "-pickle"))
    END;
    FinishedMessage("MaintainAutomaton " & Version);
  END DoIt;
  
PROCEDURE DoShowArgs(cmd: TEXT) =
  <* FATAL Thread.Alerted, Wr.Failure *>  
  BEGIN
    NL();
    Wr.PutText(arg.log, cmd);
    NL();
    NL();
  END DoShowArgs;

PROCEDURE DoLoad(rd: Rd.T; dagSize: NAT): Reduced.T =
  VAR aut: Reduced.T;
  BEGIN
    StartedMessage("loading automaton");
    aut := Reduced.Load(rd, dagSize);
    Util.PrintDoc(arg.log, aut.doc, "  |");
    FinishedMessage("loading automaton");
    RETURN aut
  END DoLoad;

PROCEDURE DoModify(
    aut: Reduced.T;
    rd: Rd.T;
    op: ModifyOp;
    noRedundant: BOOL;
    reportStep: NAT;
  ) =

  PROCEDURE NextString (* : Reduced.NextStringProc *) (
      VAR (*IO*) s: REF String; 
      VAR (*OUT*) len: NAT;
      VAR (*OUT*) add: BOOL;
    ) RAISES {Done} =
    <* FATAL Thread.Alerted, Rd.Failure *>  
    VAR c: CHAR;
    BEGIN
      TRY
        (* Read first char of line, stop if EOF: *)
        c := Rd.GetChar(rd);
        (* Decide what to do with word: *)
        CASE op OF
        | ModifyOp.Add =>
            add := TRUE;
        | ModifyOp.Remove =>
            add := FALSE;
        | ModifyOp.AddRemove =>
            add := c = '+';
            IF c # '+' AND c # '-' THEN
              Rd.UnGetChar(rd);
              BadWordMessage(arg.log, 
                "invalid add/remove line: \"" & Rd.GetLine(rd) & "\""
              );
              Process.Exit(1)
            END;
            (* Get next char, provide missing end-of-line: *)
            TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => c := '\n' END;
        ELSE
          <* ASSERT FALSE *>
        END;
        (* Put back first char of word: *)
        Rd.UnGetChar(rd);
        TRY
          encoding.ReadString(rd, s, len)
        EXCEPT
        | Done => 
            (* Should never happen: *)
            <* ASSERT FALSE *>
        | Encoding.BadChars(msg) => 
            BadWordMessage(arg.log, "badly formed word -- " & msg);
            Process.Exit(1);
        END
      EXCEPT
      | Rd.EndOfFile   =>  RAISE Done
      END
    END NextString;

  <* FATAL Abort *>
  BEGIN
    <* ASSERT aut # NIL *>
    StartedMessage(ModifyOpMessage[op]);
    aut.Build(
      next := NextString, 
      wr := arg.log,
      reportInterval := reportStep, 
      flagRedundant := NOT noRedundant
    );
    FinishedMessage(ModifyOpMessage[op])
  END DoModify;

PROCEDURE DoCheck(aut: Reduced.T; rd: Rd.T; wr: Wr.T; reportStep: NAT) =
  VAR
    c: CHAR;
    s := NEW(REF String, 100);
    len: NAT;
    badCount: NAT := 0;
    notFoundCount: NAT := 0;
    inputCount: NAT := 0;
  <* FATAL Wr.Failure, Rd.Failure *>
  BEGIN
    StartedMessage("checking words against automaton");
    TRY
      LOOP
        (* Read first char of line, stop if EOF: *)
        c := Rd.GetChar(rd);
        INC(inputCount);
        IF inputCount MOD reportStep = 0 THEN
          WrNat(arg.log, inputCount); NL();
        END;
        Rd.UnGetChar(rd);
        TRY
          encoding.ReadString(rd, s, len);
          IF NOT aut.Accepts(aut.Root(), SUBARRAY(s^, 0, len)) THEN
            INC(notFoundCount);
            FOR i := 0 TO len-1 DO
              encoding.PrintLetter(wr, s^[i])
            END;
            Wr.PutChar(wr, '\n')
          END
        EXCEPT
        | Done => 
            (* Should never happen *)
            <* ASSERT FALSE *>
        | Encoding.BadChars(msg) => 
            BadWordMessage(arg.log, "illegal word -- " & msg);
            INC(badCount);
        END;
      END
    EXCEPT
    | Rd.EndOfFile  => (* OK *)
    END;
    Wr.Flush(wr);
    NL();
    Wr.PutText(arg.log, Fmt.Int(inputCount) & " words read\n");
    Wr.PutText(arg.log, Fmt.Int(badCount) & " illegal words\n");
    Wr.PutText(arg.log, Fmt.Int(notFoundCount) & " words not accepted\n");
    FinishedMessage("checking words against automaton");
  END DoCheck;
  
PROCEDURE DoDump(aut: Reduced.T; wr: Wr.T) =
  <* FATAL Wr.Failure *> 
  BEGIN
    StartedMessage("dumping automaton");
    NL(); Util.PrintDoc(arg.log, aut.doc, "  |"); NL();
    Reduced.Dump(wr, aut);
    Wr.Close(wr);
    FinishedMessage("dumping automaton");
  END DoDump;
  
PROCEDURE MakeDoc(oldDoc, cmd, newDoc: TEXT): TEXT =
  <* FATAL Wr.Failure *> 
  BEGIN
    WITH wr = TextWr.New() DO
      Wr.PutText(wr, oldDoc);
      Wr.PutText(wr, "\n");
      Wr.PutText(wr, "--- MaintainAutomaton ");
      Wr.PutText(wr, Version);
      Wr.PutText(wr, " --- ");
      Wr.PutText(wr, Util.FmtDate(Util.GetDate()));
      Wr.PutText(wr, " ---\n");
      Wr.PutText(wr, "\n");
      Wr.PutText(wr, cmd);
      IF NOT Text.Empty(newDoc) THEN
        Wr.PutText(wr, "\n");
        Wr.PutText(wr, newDoc);
        Wr.PutText(wr, "\n");
      END;
      RETURN TextWr.ToText(wr)
    END;
  END MakeDoc;

PROCEDURE GetCommandLineArguments() =
  <* FATAL Wr.Failure *> 
  BEGIN
    WITH 
      pp = NEW(ParseParams.T).init(stderr),
      cmdWr = NEW(TextWr.T).init()
    DO
      TRY
        IF pp.keywordPresent("-help") THEN 
          pp.finish();
          Wr.PutText(stderr, HelpText);
          Wr.Flush(stderr);
          Process.Exit(0)
        ELSIF pp.keywordPresent("-version") THEN 
          pp.finish();
          Wr.PutText(stderr, Version & "\n");
          Wr.Flush(stderr);
          Process.Exit(0)
        END;
        
        Wr.PutText(cmdWr, "MaintainAutomaton");
        
        arg.log := ParamUtil.GetLogWr(pp, cmdWr);
          
        arg.dagSize := ParamUtil.GetInt(
          pp, cmdWr, "-dagSize", DefaultDagSize, 1, LAST(NAT)
        );
        
        arg.load := ParamUtil.GetFileName(pp, cmdWr, "-load");
        
        arg.modify := "";
        FOR m := FIRST(ModifyOp) TO LAST(ModifyOp) DO
          WITH f = ParamUtil.GetFileName(pp, cmdWr, ModifyOpFlag[m]) DO
            IF NOT Empty(f) THEN 
              arg.modify := f; arg.op := m; EXIT
            END
          END
        END;
        
        IF Empty(arg.load) THEN
          IF Empty(arg.modify) THEN
            Wr.PutText(stderr, "Warning: automaton will be empty")
          ELSIF arg.op = ModifyOp.Remove THEN
            Wr.PutText(stderr, "Warning: \"-remove\" not meaningful without \"-load\"")
          END
        END;

        arg.stats := ParamUtil.GetBool(pp, cmdWr, "-stats");
        
        arg.dump := ParamUtil.GetFileName(pp, cmdWr, "-dump");
        
        IF arg.dump # NIL THEN
          IF pp.keywordPresent("-comment") THEN
            arg.doc := pp.getNext();
            ParamUtil.PrintTextArg(cmdWr, "-comment", arg.doc)
          ELSE
            WITH rd = ParamUtil.GetRd(pp, cmdWr, "-commentFile") DO
              IF rd # NIL THEN 
                arg.doc := SlurpText(rd)
              ELSE
                arg.doc := ""
              END
            END
          END
        ELSE
          arg.doc := "";
        END;
        
        arg.check := ParamUtil.GetFileName(pp, cmdWr, "-check");
        arg.spell := ParamUtil.GetFileName(pp, cmdWr, "-spell");
        arg.print := ParamUtil.GetFileName(pp, cmdWr, "-print");
        arg.pickle := ParamUtil.GetFileName(pp, cmdWr, "-pickle");
        
        IF NOT Empty(arg.check) OR NOT Empty(arg.modify) THEN
          arg.reportStep := ParamUtil.GetInt(
            pp, cmdWr, "-reportStep", DefaultReportStep, 1, LAST(NAT)
          )
        ELSE
          arg.reportStep := DefaultReportStep;
        END;

        IF NOT Empty(arg.modify) THEN
          arg.noRedundant := ParamUtil.GetBool(pp, cmdWr, "-noRedundant")
        ELSE
          arg.noRedundant := FALSE
        END;

        pp.finish();
        arg.cmd := TextWr.ToText(cmdWr);
      EXCEPT
      | ParseParams.Error =>
          Wr.PutText(stderr, HelpText);
          Wr.Flush(stderr);
          Process.Exit(1)
      END
    END
  END GetCommandLineArguments;
  
PROCEDURE SlurpText(rd: Rd.T): TEXT =
  <* FATAL Rd.Failure, Wr.Failure, Rd.EndOfFile *>
  BEGIN
    WITH wr = TextWr.New() DO
      WHILE NOT Rd.EOF(rd) DO
        Wr.PutChar(wr, Rd.GetChar(rd))
      END;
      RETURN TextWr.ToText(wr)
    END;
  END SlurpText;
  
PROCEDURE StartedMessage(m: TEXT)  =
  (*
    Prints "started m" to "arg.log", with separating line. *)
  BEGIN
    NL();
    Message("started " & m);
    EVAL Util.GetExecTimes(); (* Reset the clock for FinishedMessage *)
  END StartedMessage;

PROCEDURE FinishedMessage(m: TEXT)  =
  (*
    Prints "finished m" to "arg.log", with CPU time and separating line. *)
  BEGIN
    WITH tim = Util.GetExecTimesText().total DO
      Message("finished " & m & " (" & tim & ")")
    END;
    NL();
  END FinishedMessage;

CONST MessageWidth = 76;

PROCEDURE Message(m: TEXT)  =
  (*
    Prints a separating line to "arg.log". *)
  <* FATAL Wr.Failure *>
  BEGIN
    Wr.PutText(arg.log,"===");
    IF Empty(m) THEN
      Wr.PutText(arg.log, "==")
    ELSE
      Wr.PutText(arg.log, " ");
      Wr.PutText(arg.log, m);
      Wr.PutText(arg.log, " ")
    END;
    FOR i:=1 TO MessageWidth - 5 - Text.Length(m) DO 
      Wr.PutChar(arg.log, '=')
    END;
    NL();
  END Message;

PROCEDURE NL()  =
  (* Prints a blank line to "arg.log". *)
  <* FATAL Wr.Failure *>
  BEGIN
    Wr.PutText(arg.log,"\n");
    Wr.Flush(arg.log)
  END NL;

PROCEDURE BadWordMessage(wr: Wr.T; msg: TEXT) =
  <* FATAL Wr.Failure *>
  BEGIN
    Wr.PutText(wr, "*** ");
    Wr.PutText(wr, msg);
    Wr.PutText(wr, " \n");
    Wr.Flush(wr);
  END BadWordMessage;

PROCEDURE DoWriteStats(aut: Reduced.T) =

  PROCEDURE PrintCounts()  =
    <* FATAL Wr.Failure *>
    BEGIN
      Wr.PutText(arg.log,"Counts:\n");
      Wr.PutText(arg.log,"-------\n");
      WITH counts = aut.Count(ARRAY OF Reduced.State{aut.Root()}) DO
        Reduced.PrintCounts(arg.log, counts)
      END;
      NL();
    END PrintCounts;

  PROCEDURE PrintDegreeDistr()  =
    VAR accum: NAT;
    <* FATAL Wr.Failure *>
    BEGIN
      WITH
        root = aut.Root(), 
        otrans =  NEW(REF ARRAY OF NAT, NUMBER(Letter)),
        itrans  = NEW(REF ARRAY OF NAT, root+1),
        nStates = aut.NStates(ARRAY OF Reduced.State{root}),
        rStates = FLOAT(nStates)/100.0
      DO
        FOR i := 0 TO root DO itrans^[i] := 0  END;
        FOR i := 0 TO LAST(Letter) DO otrans^[i] := 0 END;
        FOR s := Reduced.UnitState TO root DO
          IF aut.NPrefs(s)>0 THEN      
            INC(itrans^[aut.InDeg(s)]);
            INC(otrans[aut.OutDeg(s)])
          END
        END;
        accum := 0;
        Wr.PutText(arg.log,
          "\nDistribution of indegrees:" &
          "\n--------------------------\n\n" &
          "In-Trans States   %     %acc\n\n"
        );
        FOR i := 0 TO root DO
          WITH it = itrans^[i] DO
            IF it > 0 THEN
              INC(accum, it);
              Wr.PutText(arg.log,
                Fmt.Pad(Fmt.Int(i), 7) & 
                Fmt.Pad(Fmt.Int(it), 7) &
                Fmt.Pad(Fmt.Real(FLOAT(it)/rStates, Fmt.Style.Fix, 2), 7) &
                Fmt.Pad(Fmt.Real(FLOAT(accum)/rStates, Fmt.Style.Fix, 2), 7)
              );
              NL()
            END
          END
        END;
        <* ASSERT accum=nStates *>
        accum := 0;
        Wr.PutText(arg.log,
          "\n\nDistribution of outdegrees:" &
          "\n---------------------------\n\n" &
          "Out-Trans States  %     %acc\n\n"
        );
        FOR i := 0 TO LAST(Letter) DO
          WITH ot=otrans^[i] DO
            IF ot>0 THEN
              INC(accum,ot);
              Wr.PutText(arg.log,
                Fmt.Pad(Fmt.Int(i),7) & 
                Fmt.Pad(Fmt.Int(ot),7) &
                Fmt.Pad(Fmt.Real(FLOAT(ot)/rStates, Fmt.Style.Fix, 2), 7) &
                Fmt.Pad(Fmt.Real(FLOAT(accum)/rStates, Fmt.Style.Fix, 2), 7)
              );
              NL()
            END
          END
        END;
        <* ASSERT accum = nStates *>
      END;
    END PrintDegreeDistr;

  PROCEDURE PrintLengthDistr()  =
    VAR
      LenDis := ARRAY [0..MaxWordLength+1] OF NAT{0, ..};
      accum: NAT := 0;

    PROCEDURE SuffixAction(READONLY sw: String) RAISES {} =
      BEGIN
        WITH n = NUMBER(sw) DO
          IF n > MaxWordLength THEN 
            INC(LenDis[MaxWordLength+1])
          ELSE
            INC(LenDis[n])
          END
        END
      END SuffixAction;

    <* FATAL Wr.Failure, Abort *>
    BEGIN
      aut.EnumSuffs(aut.Root(), SuffixAction);
      Wr.PutText(arg.log,
        "\nDistribution of word lengths:" &
        "\n-----------------------------\n\n" &
        "Lengths   Words   %     %acc\n\n"
      );
      WITH 
        nWords = aut.NSuffs(aut.Root()),
        rWords = FLOAT(nWords)/100.0,
        hack = ARRAY BOOLEAN OF TEXT{" ", "+"}
      DO
        FOR i := 0 TO  MaxWordLength+1 DO
          WITH ldi = LenDis[i] DO
            IF ldi > 0 THEN
              INC(accum, ldi);
              Wr.PutText(arg.log,
                Fmt.Pad(Fmt.Int(i), 6) & hack[(i > MaxWordLength)] &
                Fmt.Pad(Fmt.Int(ldi), 7) &
                Fmt.Pad(Fmt.Real(FLOAT(ldi)/rWords, Fmt.Style.Fix, 2), 7) &
                Fmt.Pad(Fmt.Real(FLOAT(accum)/rWords, Fmt.Style.Fix, 2), 7)
              );
              NL()
            END
          END
        END;
        <* ASSERT accum = nWords *>
      END
    END PrintLengthDistr;

  <* FATAL Wr.Failure *>
  BEGIN
    StartedMessage("printing automaton statistics");
    PrintCounts();
    PrintDegreeDistr();
    PrintLengthDistr();
    Wr.Flush(arg.log);
    FinishedMessage("printing automaton statistics");
  END DoWriteStats;
  
PROCEDURE DoSpell(aut: Reduced.T; wr: Wr.T; reportStep: NAT) =

  VAR wordCount: NAT := 0;

  PROCEDURE EnumAction(READONLY w: String) RAISES {}=
    <* FATAL Wr.Failure, Encoding.BadString *>
    BEGIN
      encoding.PrintString(wr, w);
      Wr.PutChar(wr,'\n');
      INC(wordCount);
      IF wordCount MOD reportStep = 0  THEN
        WrNat(arg.log, wordCount);
        NL();
      END;
    END EnumAction;
  
  <* FATAL Wr.Failure, Abort *>
  BEGIN
    StartedMessage("spelling out all accepted words");
    aut.EnumSuffs(aut.Root(), EnumAction);
    Wr.Flush(wr);
    Wr.PutText(arg.log, "spelled " & Fmt.Int(wordCount) & " words\n");
    FinishedMessage("spelling out all accepted words");
  END DoSpell;
  
PROCEDURE DoPrint(aut: Reduced.T; wr: Wr.T) =
  
  <* FATAL Wr.Failure *>
  BEGIN
    StartedMessage("printing automaton");
    Util.PrintDoc(wr, aut.doc, "  |");
    Reduced.Print(wr, aut, encoding);
    Wr.Close(wr);
    FinishedMessage("printing automaton");
  END DoPrint;

PROCEDURE DoPickle(aut: Reduced.T; wr: Wr.T) =

  <* FATAL Wr.Failure, Pickle.Error *>
  BEGIN
    StartedMessage("pickling automaton");
    Pickle.Write(wr, aut);
    Wr.Close(wr);
    FinishedMessage(" pickling automaton");
  END DoPickle;

BEGIN 
  DoIt()
END MaintainAutomaton.

(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       *)
(*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         *)
(*                                                                          *)
(* This file can be freely distributed, modified, and used for any          *)
(*   non-commercial purpose, provided that this copyright and authorship    *)
(*   notice be included in any copy or derived version of this file.        *)
(*                                                                          *)
(* DISCLAIMER: This software is offered ``as is'', without any guarantee    *)
(*   as to fitness for any particular purpose.  Neither the copyright       *)
(*   holder nor the authors or their employers can be held responsible for  *)
(*   any damages that may result from its use.                              *)
(****************************************************************************)
