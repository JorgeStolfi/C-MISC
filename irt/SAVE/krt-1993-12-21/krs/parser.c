extern char *malloc(), *realloc();

# line 8 "parser.y"
#include <stdio.h>
#include "parseaux.h"

#define YYDEBUG  1
int number;
t_3d p1, p2, p3;

# line 16 "parser.y"
typedef union  {
       char *c;
       int i;
       double d;
       } YYSTYPE;
# define TINT 257
# define TFLOAT 258
# define TSTRING 259
# define TEYEP 260
# define TLOOKP 261
# define TUP 262
# define TFOV 263
# define TSCREEN 264
# define TMAXLEVEL 265
# define TLIGHT 266
# define TSURFACE 267
# define TSPHERE 268
# define TBOX 269
# define TTRIANGLE 270
# define TSUPERQ 271
# define TBACKGROUND 272
# define TOUTFILE 273
#define yyclearin yychar = -1
#define yyerrok yyerrflag = 0
extern int yychar;
extern int yyerrflag;
#ifndef YYMAXDEPTH
#define YYMAXDEPTH 150
#endif
YYSTYPE yylval, yyval;
# define YYERRCODE 256

# line 166 "parser.y"

yyerror(s)
  char *s;
  {
  fprintf( stderr,"%s\n",s);
  }

int yyexca[] ={
-1, 1,
	0, -1,
	-2, 0,
	};
# define YYNPROD 34
# define YYLAST 125
int yyact[]={

    17,    18,    19,    20,    21,    22,    24,    25,    26,    27,
    28,    29,    23,    30,    34,    33,    53,    48,    46,    45,
    44,    43,    42,    39,    38,    32,     2,    16,    31,    15,
    14,    13,    12,    11,    10,     9,     8,     7,     6,     5,
     4,     3,     1,    47,    35,    36,    37,     0,     0,    40,
    41,     0,     0,     0,     0,     0,     0,     0,    49,     0,
     0,    50,    51,    52,     0,     0,    54,    55,    56,    57,
    58,    59,    60,     0,     0,    61,    62,    63,     0,     0,
    64,    65,    66,    67,    68,    69,    70,     0,     0,     0,
     0,    71,    72,    73,    74,    75,    76,     0,    77,    78,
    79,    80,    81,    82,     0,    83,    84,    85,    86,    87,
    88,    89,    90,     0,    91,    92,    93,    94,     0,    95,
    96,    97,     0,    98,    99 };
int yypact[]={

  -260,  -260, -1000, -1000, -1000, -1000, -1000, -1000, -1000, -1000,
 -1000, -1000, -1000, -1000, -1000, -1000, -1000,  -243,  -243,  -243,
  -243,  -233,  -234,  -243,  -243,  -235,  -236,  -237,  -238,  -239,
  -242, -1000,  -243, -1000, -1000,  -243,  -243,  -243,  -241, -1000,
  -243,  -243,  -243,  -243,  -243,  -243,  -243, -1000, -1000,  -243,
  -243,  -243, -1000, -1000,  -243,  -243,  -243,  -243,  -243,  -243,
  -243, -1000, -1000, -1000, -1000,  -243,  -243,  -243,  -243,  -243,
  -243, -1000,  -243,  -243,  -243,  -243,  -243,  -243, -1000,  -243,
  -243,  -243,  -243,  -243,  -243,  -243,  -243, -1000,  -243,  -243,
  -243,  -243, -1000,  -243,  -243,  -243, -1000,  -243,  -243, -1000 };
int yypgo[]={

     0,    25,    43,    42,    26,    41,    40,    39,    38,    37,
    36,    35,    34,    33,    32,    31,    30,    29,    27 };
int yyr1[]={

     0,     3,     3,     4,     4,     4,     4,     4,     4,     4,
     4,     4,     4,     4,     4,     4,     4,     5,     6,     7,
     8,     9,    10,    11,    12,    13,    14,    15,    16,    17,
    18,     1,     1,     2 };
int yyr2[]={

     0,     4,     2,     2,     2,     2,     2,     2,     2,     2,
     2,     2,     2,     2,     2,     2,     2,     9,     9,     9,
     7,     7,     5,     9,    11,    29,    13,    17,    23,    19,
     5,     3,     3,     3 };
int yychk[]={

 -1000,    -3,    -4,    -5,    -6,    -7,    -8,    -9,   -10,   -11,
   -12,   -13,   -14,   -15,   -16,   -17,   -18,   260,   261,   262,
   263,   264,   265,   272,   266,   267,   268,   269,   270,   271,
   273,    -4,    -1,   258,   257,    -1,    -1,    -1,   257,   257,
    -1,    -1,   257,   257,   257,   257,   257,    -2,   259,    -1,
    -1,    -1,    -1,   257,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1 };
int yydef[]={

     0,    -2,     2,     3,     4,     5,     6,     7,     8,     9,
    10,    11,    12,    13,    14,    15,    16,     0,     0,     0,
     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     0,     1,     0,    31,    32,     0,     0,     0,     0,    22,
     0,     0,     0,     0,     0,     0,     0,    30,    33,     0,
     0,     0,    20,    21,     0,     0,     0,     0,     0,     0,
     0,    17,    18,    19,    23,     0,     0,     0,     0,     0,
     0,    24,     0,     0,     0,     0,     0,     0,    26,     0,
     0,     0,     0,     0,     0,     0,     0,    27,     0,     0,
     0,     0,    29,     0,     0,     0,    28,     0,     0,    25 };
typedef struct { char *t_name; int t_val; } yytoktype;
#ifndef YYDEBUG
#	define YYDEBUG	0	/* don't allow debugging */
#endif

#if YYDEBUG

yytoktype yytoks[] =
{
	"TINT",	257,
	"TFLOAT",	258,
	"TSTRING",	259,
	"TEYEP",	260,
	"TLOOKP",	261,
	"TUP",	262,
	"TFOV",	263,
	"TSCREEN",	264,
	"TMAXLEVEL",	265,
	"TLIGHT",	266,
	"TSURFACE",	267,
	"TSPHERE",	268,
	"TBOX",	269,
	"TTRIANGLE",	270,
	"TSUPERQ",	271,
	"TBACKGROUND",	272,
	"TOUTFILE",	273,
	"-unknown-",	-1	/* ends search */
};

char * yyreds[] =
{
	"-no such reduction-",
	"File : File Item",
	"File : Item",
	"Item : Eyep",
	"Item : Lookp",
	"Item : Up",
	"Item : Fov",
	"Item : Screen",
	"Item : Maxlevel",
	"Item : Background",
	"Item : Light",
	"Item : Surface",
	"Item : Sphere",
	"Item : Box",
	"Item : Triangle",
	"Item : Superq",
	"Item : Outfile",
	"Eyep : TEYEP Fnumber Fnumber Fnumber",
	"Lookp : TLOOKP Fnumber Fnumber Fnumber",
	"Up : TUP Fnumber Fnumber Fnumber",
	"Fov : TFOV Fnumber Fnumber",
	"Screen : TSCREEN TINT TINT",
	"Maxlevel : TMAXLEVEL TINT",
	"Background : TBACKGROUND Fnumber Fnumber Fnumber",
	"Light : TLIGHT Fnumber Fnumber Fnumber Fnumber",
	"Surface : TSURFACE TINT Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber",
	"Sphere : TSPHERE TINT Fnumber Fnumber Fnumber Fnumber",
	"Box : TBOX TINT Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber",
	"Triangle : TTRIANGLE TINT Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber",
	"Superq : TSUPERQ TINT Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber Fnumber",
	"Outfile : TOUTFILE String",
	"Fnumber : TFLOAT",
	"Fnumber : TINT",
	"String : TSTRING",
};
#endif /* YYDEBUG */
#line 1 "/usr/lib/yaccpar"
/*	@(#)yaccpar 1.10 89/04/04 SMI; from S5R3 1.10	*/

/*
** Skeleton parser driver for yacc output
*/

/*
** yacc user known macros and defines
*/
#define YYERROR		goto yyerrlab
#define YYACCEPT	{ free(yys); free(yyv); return(0); }
#define YYABORT		{ free(yys); free(yyv); return(1); }
#define YYBACKUP( newtoken, newvalue )\
{\
	if ( yychar >= 0 || ( yyr2[ yytmp ] >> 1 ) != 1 )\
	{\
		yyerror( "syntax error - cannot backup" );\
		goto yyerrlab;\
	}\
	yychar = newtoken;\
	yystate = *yyps;\
	yylval = newvalue;\
	goto yynewstate;\
}
#define YYRECOVERING()	(!!yyerrflag)
#ifndef YYDEBUG
#	define YYDEBUG	1	/* make debugging available */
#endif

/*
** user known globals
*/
int yydebug;			/* set to 1 to get debugging */

/*
** driver internal defines
*/
#define YYFLAG		(-1000)

/*
** static variables used by the parser
*/
static YYSTYPE *yyv;			/* value stack */
static int *yys;			/* state stack */

static YYSTYPE *yypv;			/* top of value stack */
static int *yyps;			/* top of state stack */

static int yystate;			/* current state */
static int yytmp;			/* extra var (lasts between blocks) */

int yynerrs;			/* number of errors */

int yyerrflag;			/* error recovery flag */
int yychar;			/* current input token number */


/*
** yyparse - return 0 if worked, 1 if syntax error not recovered from
*/
int
yyparse()
{
	register YYSTYPE *yypvt;	/* top of value stack for $vars */
	unsigned yymaxdepth = YYMAXDEPTH;

	/*
	** Initialize externals - yyparse may be called more than once
	*/
	yyv = (YYSTYPE*)malloc(yymaxdepth*sizeof(YYSTYPE));
	yys = (int*)malloc(yymaxdepth*sizeof(int));
	if (!yyv || !yys)
	{
		yyerror( "out of memory" );
		return(1);
	}
	yypv = &yyv[-1];
	yyps = &yys[-1];
	yystate = 0;
	yytmp = 0;
	yynerrs = 0;
	yyerrflag = 0;
	yychar = -1;

	goto yystack;
	{
		register YYSTYPE *yy_pv;	/* top of value stack */
		register int *yy_ps;		/* top of state stack */
		register int yy_state;		/* current state */
		register int  yy_n;		/* internal state number info */

		/*
		** get globals into registers.
		** branch to here only if YYBACKUP was called.
		*/
	yynewstate:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;
		goto yy_newstate;

		/*
		** get globals into registers.
		** either we just started, or we just finished a reduction
		*/
	yystack:
		yy_pv = yypv;
		yy_ps = yyps;
		yy_state = yystate;

		/*
		** top of for (;;) loop while no reductions done
		*/
	yy_stack:
		/*
		** put a state and value onto the stacks
		*/
#if YYDEBUG
		/*
		** if debugging, look up token value in list of value vs.
		** name pairs.  0 and negative (-1) are special values.
		** Note: linear search is used since time is not a real
		** consideration while debugging.
		*/
		if ( yydebug )
		{
			register int yy_i;

			(void)printf( "State %d, token ", yy_state );
			if ( yychar == 0 )
				(void)printf( "end-of-file\n" );
			else if ( yychar < 0 )
				(void)printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				(void)printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ++yy_ps >= &yys[ yymaxdepth ] )	/* room on stack? */
		{
			/*
			** reallocate and recover.  Note that pointers
			** have to be reset, or bad things will happen
			*/
			int yyps_index = (yy_ps - yys);
			int yypv_index = (yy_pv - yyv);
			int yypvt_index = (yypvt - yyv);
			yymaxdepth += YYMAXDEPTH;
			yyv = (YYSTYPE*)realloc((char*)yyv,
				yymaxdepth * sizeof(YYSTYPE));
			yys = (int*)realloc((char*)yys,
				yymaxdepth * sizeof(int));
			if (!yyv || !yys)
			{
				yyerror( "yacc stack overflow" );
				return(1);
			}
			yy_ps = yys + yyps_index;
			yy_pv = yyv + yypv_index;
			yypvt = yyv + yypvt_index;
		}
		*yy_ps = yy_state;
		*++yy_pv = yyval;

		/*
		** we have a new state - find out what to do
		*/
	yy_newstate:
		if ( ( yy_n = yypact[ yy_state ] ) <= YYFLAG )
			goto yydefault;		/* simple state */
#if YYDEBUG
		/*
		** if debugging, need to mark whether new token grabbed
		*/
		yytmp = yychar < 0;
#endif
		if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
			yychar = 0;		/* reached EOF */
#if YYDEBUG
		if ( yydebug && yytmp )
		{
			register int yy_i;

			(void)printf( "Received token " );
			if ( yychar == 0 )
				(void)printf( "end-of-file\n" );
			else if ( yychar < 0 )
				(void)printf( "-none-\n" );
			else
			{
				for ( yy_i = 0; yytoks[yy_i].t_val >= 0;
					yy_i++ )
				{
					if ( yytoks[yy_i].t_val == yychar )
						break;
				}
				(void)printf( "%s\n", yytoks[yy_i].t_name );
			}
		}
#endif /* YYDEBUG */
		if ( ( ( yy_n += yychar ) < 0 ) || ( yy_n >= YYLAST ) )
			goto yydefault;
		if ( yychk[ yy_n = yyact[ yy_n ] ] == yychar )	/*valid shift*/
		{
			yychar = -1;
			yyval = yylval;
			yy_state = yy_n;
			if ( yyerrflag > 0 )
				yyerrflag--;
			goto yy_stack;
		}

	yydefault:
		if ( ( yy_n = yydef[ yy_state ] ) == -2 )
		{
#if YYDEBUG
			yytmp = yychar < 0;
#endif
			if ( ( yychar < 0 ) && ( ( yychar = yylex() ) < 0 ) )
				yychar = 0;		/* reached EOF */
#if YYDEBUG
			if ( yydebug && yytmp )
			{
				register int yy_i;

				(void)printf( "Received token " );
				if ( yychar == 0 )
					(void)printf( "end-of-file\n" );
				else if ( yychar < 0 )
					(void)printf( "-none-\n" );
				else
				{
					for ( yy_i = 0;
						yytoks[yy_i].t_val >= 0;
						yy_i++ )
					{
						if ( yytoks[yy_i].t_val
							== yychar )
						{
							break;
						}
					}
					(void)printf( "%s\n", yytoks[yy_i].t_name );
				}
			}
#endif /* YYDEBUG */
			/*
			** look through exception table
			*/
			{
				register int *yyxi = yyexca;

				while ( ( *yyxi != -1 ) ||
					( yyxi[1] != yy_state ) )
				{
					yyxi += 2;
				}
				while ( ( *(yyxi += 2) >= 0 ) &&
					( *yyxi != yychar ) )
					;
				if ( ( yy_n = yyxi[1] ) < 0 )
					YYACCEPT;
			}
		}

		/*
		** check for syntax error
		*/
		if ( yy_n == 0 )	/* have an error */
		{
			/* no worry about speed here! */
			switch ( yyerrflag )
			{
			case 0:		/* new error */
				yyerror( "syntax error" );
				goto skip_init;
			yyerrlab:
				/*
				** get globals into registers.
				** we have a user generated syntax type error
				*/
				yy_pv = yypv;
				yy_ps = yyps;
				yy_state = yystate;
				yynerrs++;
			skip_init:
			case 1:
			case 2:		/* incompletely recovered error */
					/* try again... */
				yyerrflag = 3;
				/*
				** find state where "error" is a legal
				** shift action
				*/
				while ( yy_ps >= yys )
				{
					yy_n = yypact[ *yy_ps ] + YYERRCODE;
					if ( yy_n >= 0 && yy_n < YYLAST &&
						yychk[yyact[yy_n]] == YYERRCODE)					{
						/*
						** simulate shift of "error"
						*/
						yy_state = yyact[ yy_n ];
						goto yy_stack;
					}
					/*
					** current state has no shift on
					** "error", pop stack
					*/
#if YYDEBUG
#	define _POP_ "Error recovery pops state %d, uncovers state %d\n"
					if ( yydebug )
						(void)printf( _POP_, *yy_ps,
							yy_ps[-1] );
#	undef _POP_
#endif
					yy_ps--;
					yy_pv--;
				}
				/*
				** there is no state on stack with "error" as
				** a valid shift.  give up.
				*/
				YYABORT;
			case 3:		/* no shift yet; eat a token */
#if YYDEBUG
				/*
				** if debugging, look up token in list of
				** pairs.  0 and negative shouldn't occur,
				** but since timing doesn't matter when
				** debugging, it doesn't hurt to leave the
				** tests here.
				*/
				if ( yydebug )
				{
					register int yy_i;

					(void)printf( "Error recovery discards " );
					if ( yychar == 0 )
						(void)printf( "token end-of-file\n" );
					else if ( yychar < 0 )
						(void)printf( "token -none-\n" );
					else
					{
						for ( yy_i = 0;
							yytoks[yy_i].t_val >= 0;
							yy_i++ )
						{
							if ( yytoks[yy_i].t_val
								== yychar )
							{
								break;
							}
						}
						(void)printf( "token %s\n",
							yytoks[yy_i].t_name );
					}
				}
#endif /* YYDEBUG */
				if ( yychar == 0 )	/* reached EOF. quit */
					YYABORT;
				yychar = -1;
				goto yy_newstate;
			}
		}/* end if ( yy_n == 0 ) */
		/*
		** reduction by production yy_n
		** put stack tops, etc. so things right after switch
		*/
#if YYDEBUG
		/*
		** if debugging, print the string that is the user's
		** specification of the reduction which is just about
		** to be done.
		*/
		if ( yydebug )
			(void)printf( "Reduce by (%d) \"%s\"\n",
				yy_n, yyreds[ yy_n ] );
#endif
		yytmp = yy_n;			/* value to switch over */
		yypvt = yy_pv;			/* $vars top of value stack */
		/*
		** Look in goto table for next state
		** Sorry about using yy_state here as temporary
		** register variable, but why not, if it works...
		** If yyr2[ yy_n ] doesn't have the low order bit
		** set, then there is no action to be done for
		** this reduction.  So, no saving & unsaving of
		** registers done.  The only difference between the
		** code just after the if and the body of the if is
		** the goto yy_stack in the body.  This way the test
		** can be made before the choice of what to do is needed.
		*/
		{
			/* length of production doubled with extra bit */
			register int yy_len = yyr2[ yy_n ];

			if ( !( yy_len & 01 ) )
			{
				yy_len >>= 1;
				yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
				yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
					*( yy_ps -= yy_len ) + 1;
				if ( yy_state >= YYLAST ||
					yychk[ yy_state =
					yyact[ yy_state ] ] != -yy_n )
				{
					yy_state = yyact[ yypgo[ yy_n ] ];
				}
				goto yy_stack;
			}
			yy_len >>= 1;
			yyval = ( yy_pv -= yy_len )[1];	/* $$ = $1 */
			yy_state = yypgo[ yy_n = yyr1[ yy_n ] ] +
				*( yy_ps -= yy_len ) + 1;
			if ( yy_state >= YYLAST ||
				yychk[ yy_state = yyact[ yy_state ] ] != -yy_n )
			{
				yy_state = yyact[ yypgo[ yy_n ] ];
			}
		}
					/* save until reenter driver code */
		yystate = yy_state;
		yyps = yy_ps;
		yypv = yy_pv;
	}
	/*
	** code supplied by user is placed in this switch
	*/
	switch( yytmp )
	{
		
case 17:
# line 49 "parser.y"
{ eyep.x = yypvt[-2].d;
                  eyep.y = yypvt[-1].d;
		  eyep.z = yypvt[-0].d;
	        } break;
case 18:
# line 55 "parser.y"
{ lookp.x = yypvt[-2].d;
                  lookp.y = yypvt[-1].d;
		  lookp.z = yypvt[-0].d;
	        } break;
case 19:
# line 61 "parser.y"
{ up.x = yypvt[-2].d;
                  up.y = yypvt[-1].d;
		  up.z = yypvt[-0].d;
	        } break;
case 20:
# line 67 "parser.y"
{ hfov = yypvt[-1].d; vfov = yypvt[-0].d; } break;
case 21:
# line 70 "parser.y"
{ sizey = yypvt[-1].i; sizex = yypvt[-0].i; } break;
case 22:
# line 73 "parser.y"
{ maxlevel = yypvt[-0].i; } break;
case 23:
# line 76 "parser.y"
{ background.r = yypvt[-2].d;
                  background.g = yypvt[-1].d;
		  background.b = yypvt[-0].d;
	        } break;
case 24:
# line 83 "parser.y"
{ if (nlight == lightlim)
		      yyerror("too many lights");
		  light[nlight].bright = yypvt[-3].d;
		  light[nlight].x = yypvt[-2].d;
		  light[nlight].y = yypvt[-1].d;
		  light[nlight].z = yypvt[-0].d;
		  nlight++;
  	        } break;
case 25:
# line 97 "parser.y"
{ number = yypvt[-12].i;
		  if (number >= surfacelim)
		      yyerror("surface # too big");
		  surface[number].ar = yypvt[-11].d;
		  surface[number].ag = yypvt[-10].d;
		  surface[number].ab = yypvt[-9].d;
		  surface[number].dr = yypvt[-8].d;
		  surface[number].dg = yypvt[-7].d;
		  surface[number].db = yypvt[-6].d;
		  surface[number].sr = yypvt[-5].d;
		  surface[number].sg = yypvt[-4].d;
		  surface[number].sb = yypvt[-3].d;
		  surface[number].coef = yypvt[-2].d;
		  surface[number].refl = yypvt[-1].d;
		  surface[number].transp = yypvt[-0].d;
		  nsurface++;
	          } break;
case 26:
# line 117 "parser.y"
{ if (nobject == objectlim)
		        yyerror("too many objects");
		    maksph(yypvt[-4].i,yypvt[-3].d,yypvt[-2].d,yypvt[-1].d,yypvt[-0].d);
		    nobject++;
		  } break;
case 27:
# line 126 "parser.y"
{ if (nobject == objectlim)
		        yyerror("too many objects");
		    makbox(yypvt[-6].i,yypvt[-5].d,yypvt[-4].d,yypvt[-3].d,yypvt[-2].d,yypvt[-1].d,yypvt[-0].d);
		    nobject++;
		  } break;
case 28:
# line 136 "parser.y"
{ if (nobject == objectlim)
		        yyerror("too many objects");
		    p1.x = yypvt[-8].d; p1.y = yypvt[-7].d; p1.z = yypvt[-6].d;
		    p2.x = yypvt[-5].d; p2.y = yypvt[-4].d; p2.z = yypvt[-3].d;
		    p2.x = yypvt[-2].d; p3.y = yypvt[-1].d; p3.z = yypvt[-0].d;
		    maktri(yypvt[-9].i,&p1,&p2,&p3);
		    nobject++;
		  } break;
case 29:
# line 149 "parser.y"
{ if (nobject == objectlim)
		        yyerror("too many objects");
		    maksup(yypvt[-7].i,yypvt[-6].d,yypvt[-5].d,yypvt[-4].d,yypvt[-3].d,yypvt[-2].d,yypvt[-1].d,yypvt[-0].d);
		    nobject++;
		  } break;
case 30:
# line 156 "parser.y"
{ strcpy(outfilename, yypvt[-0].c); } break;
case 31:
# line 159 "parser.y"
{ yyval.d=yypvt[-0].d; } break;
case 32:
# line 161 "parser.y"
{ yyval.d=yypvt[-0].i; } break;
case 33:
# line 164 "parser.y"
{ yyval.c = yypvt[-0].c; } break;
	}
	goto yystack;		/* reset registers in driver code */
}
