MODULE Util;

IMPORT CPUTime, Process, Thread, Wr, Fmt;
FROM Thread IMPORT Alerted;
FROM Stdio IMPORT stderr;

PROCEDURE Message(msg: TEXT) =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    Wr.PutText(stderr, msg);
    Wr.Flush(stderr);
  END Message;

PROCEDURE Error(msg: TEXT) =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    Wr.PutText(stderr, msg);
    Wr.PutChar(stderr, '\n');
    Wr.Flush(stderr);
    Process.Exit(1)
  END Error;

PROCEDURE Digits(n: CARDINAL): CARDINAL =
  VAR d := 1; 
  BEGIN
    WHILE n >= 1000 DO n := n DIV 1000; d := d + 3 END;
    WHILE n >= 10 DO n := n DIV 10; d := d + 1 END;
    RETURN d
  END Digits;

VAR prevtime := CPUTime.Now();

PROCEDURE ResetCPUTime() =
  BEGIN
    prevtime := CPUTime.Now()
  END ResetCPUTime;

PROCEDURE WriteCPUTime(wr: Wr.T) =
  <* FATAL Wr.Failure, Alerted *>
  BEGIN
    WITH
      time = CPUTime.Now(),
      dt = ROUND(time - prevtime * 100.0D0), (* centiseconds *)
      csecs  = dt MOD 100,
      secs   = (dt DIV 100) MOD 60,
      mins   = (dt DIV 6000) MOD 60,
      hours  = (dt DIV 360000)
    DO
      Wr.PutText(wr, 
        Fmt.Pad(Fmt.Int(hours), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(mins),  2, '0') & ":" &
        Fmt.Pad(Fmt.Int(secs),  2, '0') & "." &
        Fmt.Pad(Fmt.Int(csecs), 2, '0') & "\n"
      );
      Wr.Flush(wr);
      prevtime := time;
    END
  END WriteCPUTime;

BEGIN END Util.
