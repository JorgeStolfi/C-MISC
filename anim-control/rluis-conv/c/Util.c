
#include <Util.h>



#include <CPUTime.h>
#include <Process.h>
#include <Thread.h>
#include <Wr.h>
#include <Fmt.h>
#include <Thread.h>
#include <stdio.h>

void Message(char *msg)
  <* FATAL Thread.Alerted, Wr.Failure , "??");
  {
    fprintf(stderr, "%s",  msg);
    fflush(stderr);
  } /* Message */;

void Error(char *msg)
  <* FATAL Thread.Alerted, Wr.Failure , "??");
  {
    fprintf(stderr, "%s",  msg);
    Wr.PutChar(stderr, '\n');
    fflush(stderr);
    Process.Exit(1);
  } /* Error */;

nat Digits(nat n)
  VAR d = 1; 
  {
    while (n >= 1000 ) { n = n DIV 1000; d = d + 3; }
    while (n >= 10 ) { n = n DIV 10; d = d + 1; }
    return d;
  } /* Digits */;

VAR prevtime = CPUTime.Now();

void ResetCPUTime()
  {
    prevtime = CPUTime.Now();
  } /* ResetCPUTime */;

void WriteCPUTime(FILE *wr)
  <* FATAL Wr.Failure, Alerted , "??");
  {
    { /* with*/ 
      time == CPUTime.Now(),
      dt == round(time - prevtime * 100.0D0), /* centiseconds */
      csecs == dt % 100,
      secs == (dt DIV 100) % 60,
      mins == (dt DIV 6000) % 60,
      hours == (dt DIV 360000)
    ) {
      fprintf(wr, "%s",  
        Fmt.Pad(Fmt.Int(hours), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(mins),  2, '0') & ":" &
        Fmt.Pad(Fmt.Int(secs),  2, '0') & "." &
        Fmt.Pad(Fmt.Int(csecs), 2, '0') & "\n"
      );
      fflush(wr);
      prevtime = time;
    };
  } /* WriteCPUTime */;

{; } /* Util */.
