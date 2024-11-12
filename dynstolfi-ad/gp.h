#ifndef gp_H
#define gp_H

/* gp.h - a simple 2d graphics package */
/* Luiz Henrique de Figueiredo (lhf@visgraf.impa.br) -- 16 Jan 96 */
/* Last edited on 2016-09-13 15:39:03 by stolfilocal */

#include <dv.h>

double	gpopen		(char* name);
void	gpclose		(int wait);
void	gpclear		(int wait);
void	gpflush		(void);
void	gpwait		(int t);
int	gppalette	(int c, char* name);
int	gprgb		(int c, double r, double g, double b);
int	gpcolor		(int c);
int	gpfont		(char* name);
void	gpmark		(int size, char* mark);
void	gpline		(double x1, double y1, double x2, double y2);
void	gpbox		(double xmin, double xmax, double ymin, double ymax);
void	gptri		(double x1, double y1, double x2, double y2, double x3, double y3);
void	gptext		(double x, double y, char* s, char* mode);
void	gpcircle	(double x, double y, double r);
void	gpplot		(double x, double y);
void	gpbegin		(int c);
int	gppoint		(double x, double y);
void	gpend		(void);
char*	gpevent		(int wait, double* x, double* y);
double	gpwindow	(double xmin, double xmax, double ymin, double ymax);
double	gpviewport	(double xmin, double xmax, double ymin, double ymax);
void	gpview		(double* x, double* y);
void	gpunview	(double* x, double* y);
void	gpmake		(void);

#define	gpbegin		dvbegin
#define	gpclear		dvclear
#define	gpclose		dvclose
#define	gpcolor		dvcolor
#define	gpend		dvend
#define	gpflush		dvflush
#define	gpfont		dvfont
#define	gpmark		dvmark
#define	gppalette	dvpalette
#define	gprgb		dvrgb
#define	gpwait		dvwait

#ifndef rad
#define rad(a)		((a)*(double)0.01745329252)
#endif

#ifndef round
#define round(x)	((int)((x)+(double)0.5))
#endif

#ifndef min
#define	min(x,y)	( ((x)<(y)) ? (x) : (y) )
#endif

#ifndef max
#define	max(x,y)	( ((x)>(y)) ? (x) : (y) )
#endif

#endif
