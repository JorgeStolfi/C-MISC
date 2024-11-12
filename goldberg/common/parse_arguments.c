// Hey emacs! This is -*- c++ -*-
// $Id :$

#include <unistd.h>

/////////////////////////////////////////////////////////////////////////
// function prototypes
/////////////////////////////////////////////////////////////////////////
void usage(const char *progName, const char *validOpts);

/////////////////////////////////////////////////////////////////////////
// function definitions
/////////////////////////////////////////////////////////////////////////
void parse_arguments(int argc, char *const argv[], const char *validOpts)
{
  int i;
  while ((i=getopt(argc, argv, validOpts)) != -1)
    {
      if (optarg != NULL && *optarg == '=') optarg++;
      switch (i)
	{
	case '?':
	  usage(argv[0], validOpts);
	  exit(1);
	case 'g':
	  param.gap =false;
	  break;
	case 'i':
	  param.greedy = true;
	  break;
	case 'o':
	  param.offset = atoi(optarg);
	  break;
#ifdef P5_COUNTER
	case 'p':
	  param.p5 = atoi(optarg);
	  break;
#endif
	case 'u':
	  param.globalUpdateFreq = atof(optarg);
	  if (param.globalUpdateFreq == 0) param.globalUpdateFreq = LARGE;
	  break;
      }    
    }
}

void usage(const char *progName, const char *validOpts)
{
  const char *ch;
  cerr << "Usage: " << progName;
  for(ch = validOpts; *ch != '\0'; ch++) {
    if (*ch == ':') continue;
    switch(*ch) 
      {
      case 'g':
	cerr << " -g";
	break;
      case 'i':
	cerr << " -i";
	break;
      case 'o':
	cerr << " -oX";
	break;
      case 'u':
	cerr << " -uU";
	break;
      }
  }
  cerr << "\n";
  for(ch = validOpts; *ch != '\0'; ch++) {
    if (*ch == ':') continue;
    switch(*ch) 
      {
      case 'g':
	cerr << "       -g  disables gap checking\n";
	break;
      case 'i':
	cerr << "       -i perform greedy initialization\n";
	break;
      case 'o':
	cerr << "       -oX vertex labels start at X (default "
             << param.offset << ")\n";
	break;
      case 'p':
	cerr << "       -pP (pentium only) look at pentium counter P\n";
	break;
      case 'u':
	cerr << "       -uU is global recalc parameter\n";
	cerr << "              if U=0 then global recalc is never performed\n";
	break;
      }
  }
}











