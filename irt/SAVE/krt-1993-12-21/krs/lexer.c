# include "stdio.h"
# define U(x) x
# define NLSTATE yyprevious=YYNEWLINE
# define BEGIN yybgin = yysvec + 1 +
# define INITIAL 0
# define YYLERR yysvec
# define YYSTATE (yyestate-yysvec-1)
# define YYOPTIM 1
# define YYLMAX BUFSIZ
# define output(c) putc(c,yyout)
# define input() (((yytchar=yysptr>yysbuf?U(*--yysptr):getc(yyin))==10?(yylineno++,yytchar):yytchar)==EOF?0:yytchar)
# define unput(c) {yytchar= (c);if(yytchar=='\n')yylineno--;*yysptr++=yytchar;}
# define yymore() (yymorfg=1)
# define ECHO fprintf(yyout, "%s",yytext)
# define REJECT { nstr = yyreject(); goto yyfussy;}
int yyleng; extern char yytext[];
int yymorfg;
extern char *yysptr, yysbuf[];
int yytchar;
FILE *yyin = {stdin}, *yyout = {stdout};
extern int yylineno;
struct yysvf { 
	struct yywork *yystoff;
	struct yysvf *yyother;
	int *yystops;};
struct yysvf *yyestate;
extern struct yysvf yysvec[], *yybgin;
#include <stdio.h>
#include "y.tab.h"
# define COMMENT 2
# define YYNEWLINE 10
yylex(){
int nstr; extern int yyprevious;
while((nstr = yylook()) >= 0)
yyfussy: switch(nstr){
case 0:
if(yywrap()) return(0); break;
case 1:
                     ;
break;
case 2:
                      ;
break;
case 3:
                      ;
break;
case 4:
                    {BEGIN COMMENT;}
break;
case 5:
      ;
break;
case 6:
     {if (yytext[yyleng-2] =='*')
                           BEGIN 0;}
break;
case 7:
                    {return(TEYEP);}
break;
case 8:
                   {return(TLOOKP);}
break;
case 9:
                      {return(TUP);}
break;
case 10:
                     {return(TFOV);}
break;
case 11:
                  {return(TSCREEN);}
break;
case 12:
                   {return(TLIGHT);}
break;
case 13:
                 {return(TSURFACE);}
break;
case 14:
              {return(TBACKGROUND);}
break;
case 15:
                {return(TMAXLEVEL);}
break;
case 16:
                  {return(TSPHERE);}
break;
case 17:
                     {return(TBOX);}
break;
case 18:
                {return(TTRIANGLE);}
break;
case 19:
                  {return(TSUPERQ);}
break;
case 20:
                 {return(TOUTFILE);}
break;
case 21:
                {yylval.c=yytext;
                           return(TSTRING);}
break;
case 22:
         {sscanf(yytext,"%d",&yylval.i);
                           return(TINT);}
break;
case 23:
case 24:
case 25:
{sscanf(yytext,"%F",&yylval.d);
                           return(TFLOAT);}
break;
case -1:
break;
default:
fprintf(yyout,"bad switch yylook %d",nstr);
} return(0); }
/* end of yylex */
yywrap() {return(1);}

int yyvstop[] = {
0,

2,
0,

3,
0,

1,
0,

22,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

2,
0,

3,
5,
0,

1,
0,

22,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

22,
0,

24,
0,

4,
0,

23,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

9,
21,
0,

5,
0,

6,
0,

22,
0,

24,
0,

23,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

9,
21,
0,

23,
24,
0,

25,
0,

21,
0,

17,
21,
0,

21,
0,

10,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

23,
24,
0,

25,
0,

21,
0,

17,
21,
0,

21,
0,

10,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

24,
0,

23,
0,

21,
0,

7,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

24,
0,

23,
0,

21,
0,

7,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

23,
24,
0,

21,
0,

12,
21,
0,

8,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

23,
24,
0,

21,
0,

12,
21,
0,

8,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

11,
21,
0,

16,
21,
0,

19,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

11,
21,
0,

16,
21,
0,

19,
21,
0,

21,
0,

21,
0,

21,
0,

21,
0,

20,
21,
0,

13,
21,
0,

21,
0,

21,
0,

21,
0,

20,
21,
0,

13,
21,
0,

21,
0,

21,
0,

15,
21,
0,

18,
21,
0,

21,
0,

15,
21,
0,

18,
21,
0,

21,
0,

21,
0,

14,
21,
0,

14,
21,
0,
0};
# define YYTYPE int
struct yywork { YYTYPE verify, advance; } yycrank[] = {
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	1,5,	1,6,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	1,7,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	11,45,	
1,8,	0,0,	1,9,	1,10,	
1,11,	1,12,	1,12,	1,12,	
1,12,	1,12,	1,12,	1,12,	
1,12,	1,12,	1,12,	25,63,	
27,63,	28,30,	28,63,	28,31,	
30,63,	30,65,	1,13,	1,13,	
1,13,	1,13,	1,13,	1,13,	
1,13,	1,13,	1,13,	1,13,	
1,13,	1,13,	1,13,	1,13,	
1,13,	1,13,	1,13,	1,13,	
1,13,	1,13,	1,13,	1,13,	
1,13,	1,13,	1,13,	1,13,	
29,30,	29,63,	29,64,	0,0,	
1,14,	0,0,	1,13,	1,15,	
1,13,	1,13,	1,16,	1,17,	
1,13,	1,13,	1,13,	1,13,	
1,13,	1,18,	1,19,	1,13,	
1,20,	1,13,	1,13,	1,13,	
1,21,	1,22,	1,23,	1,13,	
1,13,	1,13,	1,13,	1,13,	
3,24,	15,49,	16,51,	17,52,	
18,53,	19,55,	20,56,	22,60,	
3,25,	3,26,	18,54,	23,61,	
33,63,	33,68,	49,87,	15,50,	
50,88,	51,89,	52,90,	8,10,	
4,25,	8,12,	8,12,	8,12,	
8,12,	8,12,	8,12,	8,12,	
8,12,	8,12,	8,12,	3,27,	
9,43,	9,43,	9,43,	9,43,	
9,43,	9,43,	9,43,	9,43,	
9,43,	9,43,	3,28,	4,27,	
3,29,	3,30,	3,11,	3,31,	
10,44,	10,44,	10,44,	10,44,	
10,44,	10,44,	10,44,	10,44,	
10,44,	10,44,	4,11,	31,66,	
31,63,	31,31,	12,46,	44,82,	
3,32,	53,91,	21,57,	54,92,	
3,32,	14,48,	14,48,	14,48,	
14,48,	14,48,	14,48,	14,48,	
14,48,	14,48,	14,48,	21,58,	
32,32,	32,63,	31,67,	55,93,	
21,59,	12,47,	48,47,	56,94,	
57,95,	58,96,	59,97,	60,99,	
59,98,	83,120,	3,33,	44,82,	
87,123,	3,34,	89,124,	32,32,	
3,35,	3,36,	91,125,	32,32,	
64,66,	64,63,	64,64,	3,37,	
3,38,	4,34,	3,39,	13,13,	
4,35,	4,36,	3,40,	3,41,	
3,42,	12,47,	48,47,	4,37,	
4,38,	92,126,	4,39,	93,127,	
94,128,	83,120,	4,40,	4,41,	
4,42,	32,32,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
95,129,	96,130,	97,131,	98,132,	
13,13,	99,133,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
13,13,	13,13,	13,13,	13,13,	
24,24,	34,32,	34,63,	123,152,	
35,32,	35,63,	36,32,	36,63,	
125,153,	24,62,	126,154,	37,32,	
37,63,	127,155,	38,32,	38,63,	
103,63,	103,104,	65,63,	65,65,	
34,32,	66,63,	66,101,	35,32,	
34,32,	36,32,	128,156,	35,32,	
129,157,	36,32,	37,32,	104,63,	
104,104,	38,32,	37,32,	39,32,	
39,63,	38,32,	40,32,	40,63,	
65,100,	130,158,	24,24,	66,102,	
24,24,	24,24,	24,63,	24,24,	
131,159,	132,160,	34,32,	133,161,	
34,69,	35,32,	39,32,	36,32,	
152,174,	40,32,	39,32,	155,175,	
37,32,	40,32,	156,176,	38,32,	
24,24,	38,75,	34,70,	157,177,	
24,24,	158,178,	37,73,	36,72,	
159,179,	41,32,	41,63,	160,180,	
37,74,	134,63,	134,135,	35,71,	
42,32,	42,63,	135,63,	135,135,	
39,32,	68,63,	68,68,	40,32,	
137,63,	137,138,	161,181,	40,77,	
41,32,	67,103,	24,24,	67,103,	
41,32,	67,63,	67,104,	42,32,	
174,190,	70,32,	70,63,	42,32,	
40,78,	175,191,	39,76,	68,67,	
43,46,	40,79,	43,43,	43,43,	
43,43,	43,43,	43,43,	43,43,	
43,43,	43,43,	43,43,	43,43,	
70,32,	176,192,	41,32,	180,193,	
70,32,	69,32,	69,63,	138,63,	
138,138,	42,32,	46,83,	46,83,	
46,83,	46,83,	46,83,	46,83,	
46,83,	46,83,	46,83,	46,83,	
181,194,	41,80,	101,63,	101,101,	
69,32,	190,200,	42,81,	47,85,	
69,32,	47,85,	70,32,	46,84,	
47,86,	47,86,	47,86,	47,86,	
47,86,	47,86,	47,86,	47,86,	
47,86,	47,86,	71,32,	71,63,	
101,136,	72,32,	72,63,	73,32,	
73,63,	191,201,	74,32,	74,63,	
75,32,	75,63,	69,32,	70,106,	
162,63,	162,163,	69,105,	76,32,	
76,63,	71,32,	194,202,	46,84,	
72,32,	71,32,	73,32,	200,206,	
72,32,	74,32,	73,32,	75,32,	
100,134,	74,32,	100,134,	75,32,	
100,63,	100,135,	76,32,	77,32,	
77,63,	206,208,	76,32,	0,0,	
78,32,	78,63,	79,32,	79,63,	
102,137,	0,0,	102,137,	71,32,	
102,63,	102,138,	72,32,	0,0,	
73,32,	71,107,	77,32,	74,32,	
0,0,	75,32,	77,32,	78,32,	
73,109,	79,32,	0,0,	78,32,	
76,32,	79,32,	0,0,	80,32,	
80,63,	81,32,	81,63,	74,110,	
0,0,	72,108,	106,32,	106,63,	
108,32,	108,63,	163,63,	163,163,	
0,0,	0,0,	75,111,	0,0,	
77,32,	76,112,	80,32,	0,0,	
81,32,	78,32,	80,32,	79,32,	
81,32,	106,32,	0,0,	108,32,	
0,0,	106,32,	78,114,	108,32,	
136,162,	0,0,	136,162,	77,113,	
136,63,	136,163,	0,0,	0,0,	
79,115,	0,0,	79,116,	0,0,	
0,0,	0,0,	0,0,	0,0,	
80,32,	0,0,	81,32,	0,0,	
0,0,	0,0,	82,118,	106,32,	
82,118,	108,32,	80,117,	82,119,	
82,119,	82,119,	82,119,	82,119,	
82,119,	82,119,	82,119,	82,119,	
82,119,	84,121,	0,0,	84,121,	
0,0,	0,0,	84,122,	84,122,	
84,122,	84,122,	84,122,	84,122,	
84,122,	84,122,	84,122,	84,122,	
85,86,	85,86,	85,86,	85,86,	
85,86,	85,86,	85,86,	85,86,	
85,86,	85,86,	105,32,	105,63,	
107,32,	107,63,	0,0,	109,32,	
109,63,	0,0,	110,32,	110,63,	
111,32,	111,63,	0,0,	0,0,	
0,0,	112,32,	112,63,	0,0,	
0,0,	105,32,	0,0,	107,32,	
0,0,	105,32,	109,32,	107,32,	
0,0,	110,32,	109,32,	111,32,	
0,0,	110,32,	0,0,	111,32,	
112,32,	0,0,	0,0,	0,0,	
112,32,	113,32,	113,63,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	114,32,	114,63,	105,32,	
0,0,	107,32,	0,0,	0,0,	
109,32,	115,32,	115,63,	110,32,	
113,32,	111,32,	0,0,	105,139,	
113,32,	109,141,	112,32,	0,0,	
114,32,	0,0,	107,140,	110,142,	
114,32,	112,144,	111,143,	0,0,	
115,32,	0,0,	116,32,	116,63,	
115,32,	0,0,	0,0,	0,0,	
117,32,	117,63,	0,0,	0,0,	
0,0,	0,0,	113,32,	0,0,	
0,0,	0,0,	0,0,	0,0,	
113,145,	116,32,	114,32,	0,0,	
0,0,	116,32,	0,0,	117,32,	
114,146,	0,0,	115,32,	117,32,	
0,0,	0,0,	0,0,	0,0,	
115,147,	118,119,	118,119,	118,119,	
118,119,	118,119,	118,119,	118,119,	
118,119,	118,119,	118,119,	0,0,	
0,0,	0,0,	0,0,	116,32,	
0,0,	0,0,	0,0,	0,0,	
0,0,	117,32,	116,148,	117,149,	
120,150,	0,0,	120,150,	0,0,	
0,0,	120,151,	120,151,	120,151,	
120,151,	120,151,	120,151,	120,151,	
120,151,	120,151,	120,151,	121,122,	
121,122,	121,122,	121,122,	121,122,	
121,122,	121,122,	121,122,	121,122,	
121,122,	139,32,	139,63,	140,32,	
140,63,	0,0,	141,32,	141,63,	
142,32,	142,63,	0,0,	143,32,	
143,63,	144,32,	144,63,	0,0,	
145,32,	145,63,	146,32,	146,63,	
139,32,	0,0,	140,32,	0,0,	
139,32,	141,32,	140,32,	142,32,	
0,0,	141,32,	143,32,	142,32,	
144,32,	0,0,	143,32,	145,32,	
144,32,	146,32,	0,0,	145,32,	
0,0,	146,32,	0,0,	0,0,	
0,0,	147,32,	147,63,	0,0,	
0,0,	0,0,	139,32,	0,0,	
140,32,	0,0,	0,0,	141,32,	
0,0,	142,32,	139,164,	0,0,	
143,32,	0,0,	144,32,	0,0,	
147,32,	145,32,	143,167,	146,32,	
147,32,	148,32,	148,63,	145,169,	
144,168,	0,0,	142,166,	0,0,	
141,165,	149,32,	149,63,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	146,170,	0,0,	
148,32,	0,0,	0,0,	0,0,	
148,32,	0,0,	147,32,	0,0,	
149,32,	0,0,	164,32,	164,63,	
149,32,	150,151,	150,151,	150,151,	
150,151,	150,151,	150,151,	150,151,	
150,151,	150,151,	150,151,	165,32,	
165,63,	147,171,	166,32,	166,63,	
0,0,	164,32,	148,32,	0,0,	
148,172,	164,32,	0,0,	167,32,	
167,63,	0,0,	149,32,	0,0,	
0,0,	0,0,	165,32,	168,32,	
168,63,	166,32,	165,32,	0,0,	
0,0,	166,32,	169,32,	169,63,	
0,0,	149,173,	167,32,	0,0,	
0,0,	0,0,	167,32,	164,32,	
170,32,	170,63,	168,32,	171,32,	
171,63,	0,0,	168,32,	0,0,	
0,0,	169,32,	172,32,	172,63,	
165,32,	169,32,	0,0,	166,32,	
0,0,	0,0,	164,182,	170,32,	
0,0,	0,0,	171,32,	170,32,	
167,32,	0,0,	171,32,	173,32,	
173,63,	172,32,	0,0,	0,0,	
168,32,	172,32,	182,32,	182,63,	
0,0,	183,32,	183,63,	169,32,	
0,0,	0,0,	184,32,	184,63,	
0,0,	168,184,	173,32,	167,183,	
0,0,	170,32,	173,32,	0,0,	
171,32,	182,32,	169,185,	170,186,	
183,32,	182,32,	0,0,	172,32,	
183,32,	184,32,	0,0,	172,188,	
0,0,	184,32,	185,32,	185,63,	
0,0,	0,0,	171,187,	186,32,	
186,63,	0,0,	0,0,	0,0,	
173,32,	187,32,	187,63,	188,32,	
188,63,	0,0,	0,0,	182,32,	
173,189,	185,32,	183,32,	189,32,	
189,63,	185,32,	186,32,	184,32,	
183,196,	0,0,	186,32,	0,0,	
187,32,	184,197,	188,32,	182,195,	
187,32,	0,0,	188,32,	0,0,	
195,32,	195,63,	189,32,	0,0,	
0,0,	0,0,	189,32,	0,0,	
196,32,	196,63,	0,0,	185,32,	
0,0,	0,0,	197,32,	197,63,	
186,32,	198,32,	198,63,	195,32,	
199,32,	199,63,	187,32,	195,32,	
188,32,	203,32,	203,63,	196,32,	
204,32,	204,63,	188,198,	196,32,	
189,32,	197,32,	0,0,	0,0,	
198,32,	197,32,	0,0,	199,32,	
198,32,	205,32,	205,63,	199,32,	
203,32,	189,199,	0,0,	204,32,	
203,32,	195,32,	0,0,	204,32,	
0,0,	207,32,	207,63,	0,0,	
0,0,	196,32,	0,0,	0,0,	
205,32,	209,32,	209,63,	197,32,	
205,32,	0,0,	198,32,	0,0,	
0,0,	199,32,	196,204,	195,203,	
207,32,	0,0,	203,32,	199,205,	
207,32,	204,32,	0,0,	0,0,	
209,32,	0,0,	0,0,	0,0,	
209,32,	0,0,	0,0,	0,0,	
0,0,	203,207,	205,32,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	0,0,	0,0,	
0,0,	0,0,	207,32,	0,0,	
0,0,	0,0,	0,0,	207,209,	
0,0,	0,0,	209,32,	0,0,	
0,0};
struct yysvf yysvec[] = {
0,	0,	0,
yycrank+1,	0,		0,	
yycrank+0,	yysvec+1,	0,	
yycrank+-123,	0,		0,	
yycrank+-135,	yysvec+3,	0,	
yycrank+0,	0,		yyvstop+1,
yycrank+0,	0,		yyvstop+3,
yycrank+0,	0,		yyvstop+5,
yycrank+97,	0,		0,	
yycrank+108,	yysvec+8,	0,	
yycrank+124,	0,		0,	
yycrank+1,	0,		0,	
yycrank+140,	yysvec+8,	yyvstop+7,
yycrank+189,	0,		yyvstop+9,
yycrank+145,	0,		0,	
yycrank+28,	yysvec+13,	yyvstop+11,
yycrank+5,	yysvec+13,	yyvstop+13,
yycrank+16,	yysvec+13,	yyvstop+15,
yycrank+23,	yysvec+13,	yyvstop+17,
yycrank+32,	yysvec+13,	yyvstop+19,
yycrank+13,	yysvec+13,	yyvstop+21,
yycrank+91,	yysvec+13,	yyvstop+23,
yycrank+17,	yysvec+13,	yyvstop+25,
yycrank+23,	yysvec+13,	yyvstop+27,
yycrank+-311,	0,		0,	
yycrank+-12,	yysvec+24,	yyvstop+29,
yycrank+0,	0,		yyvstop+31,
yycrank+-13,	yysvec+24,	yyvstop+34,
yycrank+-15,	yysvec+24,	0,	
yycrank+-46,	yysvec+24,	0,	
yycrank+-17,	yysvec+24,	0,	
yycrank+-137,	yysvec+24,	yyvstop+36,
yycrank+-158,	yysvec+24,	yyvstop+38,
yycrank+-89,	yysvec+24,	0,	
yycrank+-267,	yysvec+24,	yyvstop+40,
yycrank+-270,	yysvec+24,	yyvstop+42,
yycrank+-272,	yysvec+24,	yyvstop+44,
yycrank+-277,	yysvec+24,	yyvstop+46,
yycrank+-280,	yysvec+24,	yyvstop+48,
yycrank+-301,	yysvec+24,	yyvstop+50,
yycrank+-304,	yysvec+24,	yyvstop+52,
yycrank+-339,	yysvec+24,	yyvstop+54,
yycrank+-346,	yysvec+24,	yyvstop+56,
yycrank+374,	0,		yyvstop+58,
yycrank+118,	yysvec+10,	yyvstop+60,
yycrank+0,	0,		yyvstop+62,
yycrank+394,	0,		yyvstop+64,
yycrank+416,	0,		0,	
yycrank+141,	yysvec+14,	0,	
yycrank+39,	yysvec+13,	yyvstop+66,
yycrank+20,	yysvec+13,	yyvstop+68,
yycrank+40,	yysvec+13,	yyvstop+70,
yycrank+24,	yysvec+13,	yyvstop+72,
yycrank+86,	yysvec+13,	yyvstop+74,
yycrank+80,	yysvec+13,	yyvstop+76,
yycrank+87,	yysvec+13,	yyvstop+78,
yycrank+95,	yysvec+13,	yyvstop+80,
yycrank+98,	yysvec+13,	yyvstop+82,
yycrank+109,	yysvec+13,	yyvstop+84,
yycrank+102,	yysvec+13,	yyvstop+86,
yycrank+110,	yysvec+13,	yyvstop+88,
yycrank+0,	yysvec+13,	yyvstop+90,
yycrank+0,	0,		yyvstop+93,
yycrank+0,	0,		yyvstop+95,
yycrank+-182,	yysvec+24,	yyvstop+97,
yycrank+-283,	yysvec+24,	yyvstop+99,
yycrank+-286,	yysvec+24,	yyvstop+101,
yycrank+-362,	yysvec+24,	0,	
yycrank+-350,	yysvec+24,	0,	
yycrank+-391,	yysvec+24,	yyvstop+103,
yycrank+-367,	yysvec+24,	yyvstop+105,
yycrank+-428,	yysvec+24,	yyvstop+107,
yycrank+-431,	yysvec+24,	yyvstop+109,
yycrank+-433,	yysvec+24,	yyvstop+111,
yycrank+-436,	yysvec+24,	yyvstop+113,
yycrank+-438,	yysvec+24,	yyvstop+115,
yycrank+-445,	yysvec+24,	yyvstop+117,
yycrank+-465,	yysvec+24,	yyvstop+119,
yycrank+-470,	yysvec+24,	yyvstop+121,
yycrank+-472,	yysvec+24,	yyvstop+123,
yycrank+-497,	yysvec+24,	yyvstop+125,
yycrank+-499,	yysvec+24,	yyvstop+127,
yycrank+555,	0,		0,	
yycrank+148,	yysvec+46,	yyvstop+130,
yycrank+570,	0,		0,	
yycrank+580,	0,		0,	
yycrank+0,	yysvec+85,	yyvstop+133,
yycrank+113,	yysvec+13,	yyvstop+135,
yycrank+0,	yysvec+13,	yyvstop+137,
yycrank+110,	yysvec+13,	yyvstop+140,
yycrank+0,	yysvec+13,	yyvstop+142,
yycrank+122,	yysvec+13,	yyvstop+145,
yycrank+138,	yysvec+13,	yyvstop+147,
yycrank+139,	yysvec+13,	yyvstop+149,
yycrank+146,	yysvec+13,	yyvstop+151,
yycrank+179,	yysvec+13,	yyvstop+153,
yycrank+180,	yysvec+13,	yyvstop+155,
yycrank+181,	yysvec+13,	yyvstop+157,
yycrank+181,	yysvec+13,	yyvstop+159,
yycrank+188,	yysvec+13,	yyvstop+161,
yycrank+-461,	yysvec+24,	0,	
yycrank+-407,	yysvec+24,	yyvstop+163,
yycrank+-477,	yysvec+24,	0,	
yycrank+-281,	yysvec+24,	0,	
yycrank+-296,	yysvec+24,	yyvstop+166,
yycrank+-592,	yysvec+24,	yyvstop+168,
yycrank+-504,	yysvec+24,	yyvstop+170,
yycrank+-594,	yysvec+24,	yyvstop+173,
yycrank+-506,	yysvec+24,	yyvstop+175,
yycrank+-597,	yysvec+24,	yyvstop+178,
yycrank+-600,	yysvec+24,	yyvstop+180,
yycrank+-602,	yysvec+24,	yyvstop+182,
yycrank+-607,	yysvec+24,	yyvstop+184,
yycrank+-631,	yysvec+24,	yyvstop+186,
yycrank+-639,	yysvec+24,	yyvstop+188,
yycrank+-647,	yysvec+24,	yyvstop+190,
yycrank+-668,	yysvec+24,	yyvstop+192,
yycrank+-674,	yysvec+24,	yyvstop+194,
yycrank+701,	0,		0,	
yycrank+0,	yysvec+118,	yyvstop+196,
yycrank+729,	0,		0,	
yycrank+739,	0,		0,	
yycrank+0,	yysvec+121,	yyvstop+198,
yycrank+212,	yysvec+13,	yyvstop+200,
yycrank+0,	yysvec+13,	yyvstop+202,
yycrank+204,	yysvec+13,	yyvstop+205,
yycrank+210,	yysvec+13,	yyvstop+207,
yycrank+224,	yysvec+13,	yyvstop+209,
yycrank+233,	yysvec+13,	yyvstop+211,
yycrank+239,	yysvec+13,	yyvstop+213,
yycrank+239,	yysvec+13,	yyvstop+215,
yycrank+246,	yysvec+13,	yyvstop+217,
yycrank+264,	yysvec+13,	yyvstop+219,
yycrank+253,	yysvec+13,	yyvstop+221,
yycrank+-342,	yysvec+24,	0,	
yycrank+-347,	yysvec+24,	yyvstop+223,
yycrank+-533,	yysvec+24,	0,	
yycrank+-353,	yysvec+24,	0,	
yycrank+-392,	yysvec+24,	yyvstop+225,
yycrank+-751,	yysvec+24,	yyvstop+227,
yycrank+-753,	yysvec+24,	yyvstop+229,
yycrank+-756,	yysvec+24,	yyvstop+232,
yycrank+-758,	yysvec+24,	yyvstop+234,
yycrank+-761,	yysvec+24,	yyvstop+236,
yycrank+-763,	yysvec+24,	yyvstop+238,
yycrank+-766,	yysvec+24,	yyvstop+240,
yycrank+-768,	yysvec+24,	yyvstop+242,
yycrank+-795,	yysvec+24,	yyvstop+244,
yycrank+-819,	yysvec+24,	yyvstop+246,
yycrank+-827,	yysvec+24,	yyvstop+248,
yycrank+849,	0,		0,	
yycrank+0,	yysvec+150,	yyvstop+250,
yycrank+254,	yysvec+13,	yyvstop+253,
yycrank+0,	yysvec+13,	yyvstop+255,
yycrank+0,	yysvec+13,	yyvstop+258,
yycrank+253,	yysvec+13,	yyvstop+261,
yycrank+266,	yysvec+13,	yyvstop+263,
yycrank+269,	yysvec+13,	yyvstop+265,
yycrank+280,	yysvec+13,	yyvstop+267,
yycrank+271,	yysvec+13,	yyvstop+269,
yycrank+288,	yysvec+13,	yyvstop+271,
yycrank+299,	yysvec+13,	yyvstop+273,
yycrank+-441,	yysvec+24,	0,	
yycrank+-507,	yysvec+24,	yyvstop+275,
yycrank+-848,	yysvec+24,	yyvstop+278,
yycrank+-861,	yysvec+24,	yyvstop+280,
yycrank+-864,	yysvec+24,	yyvstop+283,
yycrank+-873,	yysvec+24,	yyvstop+286,
yycrank+-881,	yysvec+24,	yyvstop+288,
yycrank+-888,	yysvec+24,	yyvstop+290,
yycrank+-898,	yysvec+24,	yyvstop+292,
yycrank+-901,	yysvec+24,	yyvstop+294,
yycrank+-908,	yysvec+24,	yyvstop+296,
yycrank+-925,	yysvec+24,	yyvstop+298,
yycrank+301,	yysvec+13,	yyvstop+300,
yycrank+316,	yysvec+13,	yyvstop+302,
yycrank+332,	yysvec+13,	yyvstop+304,
yycrank+0,	yysvec+13,	yyvstop+306,
yycrank+0,	yysvec+13,	yyvstop+309,
yycrank+0,	yysvec+13,	yyvstop+312,
yycrank+334,	yysvec+13,	yyvstop+315,
yycrank+344,	yysvec+13,	yyvstop+317,
yycrank+-932,	yysvec+24,	yyvstop+319,
yycrank+-935,	yysvec+24,	yyvstop+321,
yycrank+-940,	yysvec+24,	yyvstop+323,
yycrank+-964,	yysvec+24,	yyvstop+325,
yycrank+-969,	yysvec+24,	yyvstop+328,
yycrank+-975,	yysvec+24,	yyvstop+331,
yycrank+-977,	yysvec+24,	yyvstop+334,
yycrank+-985,	yysvec+24,	yyvstop+336,
yycrank+340,	yysvec+13,	yyvstop+338,
yycrank+373,	yysvec+13,	yyvstop+340,
yycrank+0,	yysvec+13,	yyvstop+342,
yycrank+0,	yysvec+13,	yyvstop+345,
yycrank+393,	yysvec+13,	yyvstop+348,
yycrank+-1002,	yysvec+24,	yyvstop+350,
yycrank+-1010,	yysvec+24,	yyvstop+352,
yycrank+-1016,	yysvec+24,	yyvstop+354,
yycrank+-1019,	yysvec+24,	yyvstop+357,
yycrank+-1022,	yysvec+24,	yyvstop+360,
yycrank+389,	yysvec+13,	yyvstop+362,
yycrank+0,	yysvec+13,	yyvstop+364,
yycrank+0,	yysvec+13,	yyvstop+367,
yycrank+-1027,	yysvec+24,	yyvstop+370,
yycrank+-1030,	yysvec+24,	yyvstop+372,
yycrank+-1043,	yysvec+24,	yyvstop+375,
yycrank+413,	yysvec+13,	yyvstop+378,
yycrank+-1055,	yysvec+24,	yyvstop+380,
yycrank+0,	yysvec+13,	yyvstop+382,
yycrank+-1063,	yysvec+24,	yyvstop+385,
0,	0,	0};
struct yywork *yytop = yycrank+1158;
struct yysvf *yybgin = yysvec+1;
char yymatch[] = {
00  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,012 ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,01  ,01  ,'+' ,01  ,'-' ,'.' ,012 ,
'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,'0' ,
'0' ,'0' ,01  ,01  ,01  ,01  ,01  ,01  ,
01  ,'A' ,'A' ,'A' ,'A' ,'E' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,01  ,01  ,01  ,01  ,'_' ,
01  ,'A' ,'A' ,'A' ,'A' ,'E' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,'A' ,
'A' ,'A' ,'A' ,01  ,01  ,01  ,01  ,01  ,
0};
char yyextra[] = {
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,
0};
#ifndef lint
static	char ncform_sccsid[] = "@(#)ncform 1.6 88/02/08 SMI"; /* from S5R2 1.2 */
#endif

int yylineno =1;
# define YYU(x) x
# define NLSTATE yyprevious=YYNEWLINE
char yytext[YYLMAX];
struct yysvf *yylstate [YYLMAX], **yylsp, **yyolsp;
char yysbuf[YYLMAX];
char *yysptr = yysbuf;
int *yyfnd;
extern struct yysvf *yyestate;
int yyprevious = YYNEWLINE;
yylook(){
	register struct yysvf *yystate, **lsp;
	register struct yywork *yyt;
	struct yysvf *yyz;
	int yych, yyfirst;
	struct yywork *yyr;
# ifdef LEXDEBUG
	int debug;
# endif
	char *yylastch;
	/* start off machines */
# ifdef LEXDEBUG
	debug = 0;
# endif
	yyfirst=1;
	if (!yymorfg)
		yylastch = yytext;
	else {
		yymorfg=0;
		yylastch = yytext+yyleng;
		}
	for(;;){
		lsp = yylstate;
		yyestate = yystate = yybgin;
		if (yyprevious==YYNEWLINE) yystate++;
		for (;;){
# ifdef LEXDEBUG
			if(debug)fprintf(yyout,"state %d\n",yystate-yysvec-1);
# endif
			yyt = yystate->yystoff;
			if(yyt == yycrank && !yyfirst){  /* may not be any transitions */
				yyz = yystate->yyother;
				if(yyz == 0)break;
				if(yyz->yystoff == yycrank)break;
				}
			*yylastch++ = yych = input();
			yyfirst=0;
		tryagain:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"char ");
				allprint(yych);
				putchar('\n');
				}
# endif
			yyr = yyt;
			if ( (int)yyt > (int)yycrank){
				yyt = yyr + yych;
				if (yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				}
# ifdef YYOPTIM
			else if((int)yyt < (int)yycrank) {		/* r < yycrank */
				yyt = yyr = yycrank+(yycrank-yyt);
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"compressed state\n");
# endif
				yyt = yyt + yych;
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transitions */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				yyt = yyr + YYU(yymatch[yych]);
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"try fall back character ");
					allprint(YYU(yymatch[yych]));
					putchar('\n');
					}
# endif
				if(yyt <= yytop && yyt->verify+yysvec == yystate){
					if(yyt->advance+yysvec == YYLERR)	/* error transition */
						{unput(*--yylastch);break;}
					*lsp++ = yystate = yyt->advance+yysvec;
					goto contin;
					}
				}
			if ((yystate = yystate->yyother) && (yyt= yystate->yystoff) != yycrank){
# ifdef LEXDEBUG
				if(debug)fprintf(yyout,"fall back to state %d\n",yystate-yysvec-1);
# endif
				goto tryagain;
				}
# endif
			else
				{unput(*--yylastch);break;}
		contin:
# ifdef LEXDEBUG
			if(debug){
				fprintf(yyout,"state %d char ",yystate-yysvec-1);
				allprint(yych);
				putchar('\n');
				}
# endif
			;
			}
# ifdef LEXDEBUG
		if(debug){
			fprintf(yyout,"stopped at %d with ",*(lsp-1)-yysvec-1);
			allprint(yych);
			putchar('\n');
			}
# endif
		while (lsp-- > yylstate){
			*yylastch-- = 0;
			if (*lsp != 0 && (yyfnd= (*lsp)->yystops) && *yyfnd > 0){
				yyolsp = lsp;
				if(yyextra[*yyfnd]){		/* must backup */
					while(yyback((*lsp)->yystops,-*yyfnd) != 1 && lsp > yylstate){
						lsp--;
						unput(*yylastch--);
						}
					}
				yyprevious = YYU(*yylastch);
				yylsp = lsp;
				yyleng = yylastch-yytext+1;
				yytext[yyleng] = 0;
# ifdef LEXDEBUG
				if(debug){
					fprintf(yyout,"\nmatch ");
					sprint(yytext);
					fprintf(yyout," action %d\n",*yyfnd);
					}
# endif
				return(*yyfnd++);
				}
			unput(*yylastch);
			}
		if (yytext[0] == 0  /* && feof(yyin) */)
			{
			yysptr=yysbuf;
			return(0);
			}
		yyprevious = yytext[0] = input();
		if (yyprevious>0)
			output(yyprevious);
		yylastch=yytext;
# ifdef LEXDEBUG
		if(debug)putchar('\n');
# endif
		}
	}
yyback(p, m)
	int *p;
{
if (p==0) return(0);
while (*p)
	{
	if (*p++ == m)
		return(1);
	}
return(0);
}
	/* the following are only used in the lex library */
yyinput(){
	return(input());
	}
yyoutput(c)
  int c; {
	output(c);
	}
yyunput(c)
   int c; {
	unput(c);
	}
