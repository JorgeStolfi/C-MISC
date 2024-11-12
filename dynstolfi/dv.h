#ifndef dv_H
#define dv_H

/* dv.h - graphics primitives for {gp.h} */
/* Luiz Henrique de Figueiredo (lhf@visgraf.impa.br) -- 16 Jan 96 */
/* Last edited on 2008-04-24 15:28:25 by stolfi */

#define real double

typedef struct Box
{
 real xmin;
 real xmax;
 real ymin;
 real ymax;
 real xu;
 real yu;
} Box;

Box*	dvopen		(char* name);
void    dvbusy          (int busy);
void	dvclose		(int wait);
void	dvclear		(int wait);
void	dvflush		(void);
void	dvwait		(int t);
int	dvpalette	(int c, char* name);
int	dvrgb		(int c, real r, real g, real b);
int	dvcolor		(int c);
int	dvfont		(char* name);
void	dvmark		(int size, char* mark);
void	dvclip		(int xmin, int xmax, int ymin, int ymax);
void	dvline		(int x1, int y1, int x2, int y2);
void	dvbbox		(int xmin, int xmax, int ymin, int ymax);
void	dvbox		(int xmin, int xmax, int ymin, int ymax);
void	dvtri		(int x1, int y1, int x2, int y2, int x3, int y3);
void	dvtext		(int x, int y, char* s, char* mode);
void	dvcircle	(int x, int y, int r);
void	dvplot		(int x, int y);
void	dvbegin		(int c);
int	dvpoint		(int x, int y);
void	dvend		(void);
char*	dvevent		(int wait, int* x, int* y);
void    dvdump          (int n);

void dvdoublebuffer(int on);
void dvswapbuffers(void);
void dvfrontbuffer(void);
void dvbackbuffer(void);
void dvbufferarea(int xmin, int xmax, int ymin, int ymax);

#endif
