/* See SOeps.h */
/* Last edited on 2017-06-21 05:32:54 by stolfilocal */

#include <SOeps.h>

#include <affirm.h>
#include <nat.h>

#include <stdio.h>
#include <math.h>

/*** PROTOTYPES FOR INTERNAL FUNCTIONS ***/

void eps_aux_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double r, double g, double b,
    char *operator
  );
  
void eps_aux_polygon(
    FILE *psfile,
    double x[], double y[],
    int npoints,
    double r, double g, double b,
    char *operator
  );

void eps_aux_circle(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b,
    char *operator
  );
  
void eps_aux_dot(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b,
    char *operator
  );
void eps_aux_lune(
    FILE *psfile,
    double xc, double yc, double radius, double tilt,
    double r, double g, double b,
    char *operator
  );
  
void eps_aux_slice(
    FILE *psfile,
    double xc, double yc, double radius, double start, double stop,
    double r, double g, double b,
    char *operator
  );
  
void eps_define_procs(FILE *psfile);
void eps_define_caption_procs(FILE *psfile);

void eps_save_scales(
    double xmin, double xmax,
    double ymin, double ymax,
    double hmin, double hmax,
    double vmin, double vmax
  );
  
void eps_setup_page_state(FILE *psfile, int xn, int yn);
void eps_setup_caption_data(FILE *psfile);
void eps_aux_end_page(FILE *psfile);

void eps_put_text(FILE *psfile, char *text, char *newline);

/*** IMPLEMENTATIONS ***/

static double ps_hmin, ps_hmax;
static double ps_vmin, ps_vmax;

static double ps_xmin, ps_xmax, ps_xscale;
static double ps_ymin, ps_ymax, ps_yscale;

void eps_begin_figure(
    FILE *psfile,
    double xmin, double xmax,
    double ymin, double ymax,

    double hmin, double hmax,
    double vmin, double vmax,

    int xn, int yn
  )
  {
    char *date = today();

    fprintf(stderr, "date = %s\n", date);
    
    eps_save_scales(
      xmin, xmax, ymin, ymax,
      hmin, hmax, vmin, vmax
    );
    
    fprintf(psfile, "%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(psfile, "%%%%CreationDate: %s\n", date);
    fprintf(psfile, "%%%%BoundingBox: %.0f %.0f %.0f %.0f\n", 
      floor(ps_hmin), floor(ps_vmin), ceil(ps_hmax), ceil(ps_vmax)
    );
    fprintf(psfile, "%%%%Pages: 1\n");
    fprintf(psfile, "%%%%EndComments\n");
    
    fprintf(psfile, "/$maindict 6400 dict def \n");
    fprintf(psfile, "$maindict begin\n");
    
    eps_define_procs(psfile);
    eps_define_caption_procs(psfile);
    
    fprintf(psfile, "end\n");
    fprintf(psfile, "%%%%EndProlog\n");

    fprintf(psfile, "%%%%Page: 1\n");
    fprintf(psfile, "$maindict begin\n");
    fprintf(psfile, "/savedstate save def\n");

    eps_setup_page_state(psfile, xn, yn);
    eps_setup_caption_data(psfile);
    
    fflush(psfile);
  }
  
void eps_end_figure(FILE *psfile)
  {
    eps_aux_end_page(psfile);
    fprintf(psfile, "%%%%Trailer\n");
    fflush(psfile);
  }

void eps_begin_document(FILE *psfile)
  {
    fprintf(psfile, "%%!PS-Adobe-2.0\n");
    fprintf(psfile, "%%%%Pages: (atend)\n");
    fprintf(psfile, "%%%%EndComments\n");
    
    fprintf(psfile, "/$maindict 6400 dict def \n");
    fprintf(psfile, "$maindict begin\n");
    
    eps_define_procs(psfile);
    eps_define_caption_procs(psfile);
    
    fprintf(psfile, "end\n");
    fprintf(psfile, "%%%%EndProlog\n");

    fflush(psfile);
  }

void eps_end_document(FILE *psfile, int npages)
  {
    fprintf(psfile, "%%%%Trailer\n");
    fprintf(psfile, "%%%%Pages: %d\n", npages);
    fflush(psfile);
  }

void eps_begin_page(
    FILE *psfile,
    int page,
    double xmin, double xmax,
    double ymin, double ymax,

    double hmin, double hmax,
    double vmin, double vmax,

    int xn, int yn
  )
  {
    char *date = today();

    fprintf(stderr, "date = %s\n", date);
    
    eps_save_scales(
      xmin, xmax, ymin, ymax,
      hmin, hmax, vmin, vmax
    );

    fprintf(psfile, "%%%%Page: %d %d\n", page, page);
    fprintf(psfile, "\n");

    fprintf(psfile, "$maindict begin\n");
    fprintf(psfile, "/savedstate save def\n");

    fprintf(psfile, "%% Print date:\n");
    fprintf(psfile, "gsave\n");
    fprintf(psfile, "  /Courier findfont\n");
    fprintf(psfile, "  10 scalefont setfont\n");
    fprintf(psfile, "  80 18 moveto\n");
    fprintf(psfile, "  (%s   page %d) show\n", date, page);
    fprintf(psfile, "grestore\n");
    fprintf(psfile, "\n");

    eps_setup_page_state(psfile, xn, yn);
    eps_setup_caption_data(psfile);
  }

void eps_end_page(FILE *psfile)
  {
    fprintf(psfile, "showpage\n");
    fprintf(psfile, "%%%%EndPage\n");
    fprintf(psfile, "\n");
    eps_aux_end_page(psfile);
  }

void eps_save_scales(
    double xmin, double xmax,
    double ymin, double ymax,

    double hmin, double hmax,
    double vmin, double vmax
  )
  {
    ps_hmin = hmin;
    ps_hmax = hmax;
    ps_xmin = xmin;
    ps_xmax = xmax;
    ps_xscale = (hmax - hmin)/(xmax - xmin);

    ps_vmin = vmin;
    ps_vmax = vmax;
    ps_ymin = ymin;
    ps_ymax = ymax;
    ps_yscale = (vmax - vmin)/(ymax - ymin);

    affirm (
      (fabs(ps_xscale - ps_yscale)/fabs(ps_xscale + ps_yscale) < 0.0001),
      "eps_begin_figure: unequal Postscript scales"
    );
  }

void eps_define_procs(FILE *psfile)
  {
    fprintf(psfile, "%% Rectangle draw operator:\n");
    fprintf(psfile, "%%   /xlo/ /xhi/ /ylo/ /yhi/ recd --> \n");
    fprintf(psfile, "/recd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "    3 index 2 index moveto\n");
    fprintf(psfile, "    2 index 2 index lineto\n");
    fprintf(psfile, "    2 index 1 index lineto\n");
    fprintf(psfile, "    3 index 1 index lineto\n");
    fprintf(psfile, "    pop pop pop pop\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Rectangle fill operator:\n");
    fprintf(psfile, "%%   /xlo/ /xhi/ /ylo/ /yhi/ /r/ /g/ /b/ recf --> \n");
    fprintf(psfile, "/recf\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    setrgbcolor\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    3 index 2 index moveto\n");
    fprintf(psfile, "    2 index 2 index lineto\n");
    fprintf(psfile, "    2 index 1 index lineto\n");
    fprintf(psfile, "    3 index 1 index lineto\n");
    fprintf(psfile, "    pop pop pop pop\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    fill\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Rectangle fill and stroke operator:\n");
    fprintf(psfile, "%%   /xlo/ /xhi/ /ylo/ /yhi/ /r/ /g/ /b/ recfd --> \n");
    fprintf(psfile, "/recfd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    6 index 5 index moveto\n");
    fprintf(psfile, "    5 index 5 index lineto\n");
    fprintf(psfile, "    5 index 4 index lineto\n");
    fprintf(psfile, "    6 index 4 index lineto\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    gsave setrgbcolor fill grestore\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "    pop pop pop pop\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Polygon draw operator:\n");
    fprintf(psfile, "%%   /x1/ /y1/ .. /xn/ /yn/ /n/ pold --> \n");
    fprintf(psfile, "/pold\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    3 1 roll moveto\n");
    fprintf(psfile, "    1 sub 1 1  3 2 roll { pop lineto } for\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Polygon fill operator:\n");
    fprintf(psfile, "%%   /x1/ /y1/ .. /xn/ /yn/ /n/ /r/ /g/ /b/ polf --> \n");
    fprintf(psfile, "/polf\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    setrgbcolor\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    3 1 roll moveto\n");
    fprintf(psfile, "    1 sub 1 1  3 2 roll { pop lineto } for\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    fill\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Polygon fill and draw operator:\n");
    fprintf(psfile, "%%   /x1/ /y1/ .. /xn/ /yn/ /n/ /r/ /g/ /b/ polfd --> \n");
    fprintf(psfile, "/polfd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    6 4 roll moveto\n");
    fprintf(psfile, "    4 3 roll 1 sub\n");
    fprintf(psfile, "    %% -- /x1/ /y1/ .. /x[n-1]/ /y[n-1]/ /r/ /g/ /b/ /n-1/\n");
    fprintf(psfile, "    1 1  3 2 roll { pop 5 3 roll lineto } for\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    gsave setrgbcolor fill grestore\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Circle fill operator:\n");
    fprintf(psfile, "%%   /x/ /y/ /radius/ /r/ /g/ /b/ cirf --> \n");
    fprintf(psfile, "/cirf\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    setrgbcolor\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "      0 360 arc\n");
    fprintf(psfile, "      closepath\n");
    fprintf(psfile, "    fill\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Circle draw operator:\n");
    fprintf(psfile, "%%   /x/ /y/ /radius/ cird --> \n");
    fprintf(psfile, "/cird\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "      0 360 arc\n");
    fprintf(psfile, "      closepath\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Circle fill and draw operator:\n");
    fprintf(psfile, "%%   /x/ /y/ /radius/ /r/ /g/ /b/ cirfd --> \n");
    fprintf(psfile, "/cirfd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    6 3 roll\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "      0 360 arc\n");
    fprintf(psfile, "      closepath\n");
    fprintf(psfile, "      gsave setrgbcolor fill grestore\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Lune fill and draw operator:\n");
    fprintf(psfile, "%%   /xc/ /yc/ /rad/ /tilt/ /r/ /g/ /b/ lunefd --> \n");
    fprintf(psfile, "/lunefd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    7 3 roll\n");
    fprintf(psfile, "    %% --- r, g, b, xc, yc, rad, tilt\n");
    fprintf(psfile, "    4 2 roll\n");
    fprintf(psfile, "    %% --- r, g, b, rad, tilt, xc, yc\n");
    fprintf(psfile, "    translate\n");
    fprintf(psfile, "    %% --- r, g, b, rad, tilt\n");
    fprintf(psfile, "    rotate\n");
    fprintf(psfile, "    %% --- r, g, b, rad \n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "      dup\n");
    fprintf(psfile, "      %% --- r, g, b, rad, rad\n");
    fprintf(psfile, "      dup neg 0\n");
    fprintf(psfile, "      %% --- r, g, b, rad, rad, -rad, 0\n");
    fprintf(psfile, "      3 2 roll 2 mul\n");
    fprintf(psfile, "      %% --- r, g, b, rad, -rad, 0, 2*rad\n");
    fprintf(psfile, "      -60 60 arc\n");
    fprintf(psfile, "      %% --- r, g, b, rad\n");
    fprintf(psfile, "      dup 0\n");
    fprintf(psfile, "      %% --- r, g, b, rad, rad, 0\n");
    fprintf(psfile, "      3 2 roll 2 mul\n");
    fprintf(psfile, "      %% --- r, g, b, rad, 0, 2*rad\n");
    fprintf(psfile, "      120 240 arc\n");
    fprintf(psfile, "      closepath\n");
    fprintf(psfile, "      gsave setrgbcolor fill grestore\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Slice fill and draw operator:\n");
    fprintf(psfile, "%%   /xc/ /yc/ /rad/ /start/ /stop/ /r/ /g/ /b/ slicefd --> \n");
    fprintf(psfile, "/slicefd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    8 3 roll\n");
    fprintf(psfile, "    %% --- r, g, b, xc, yc, rad, start, stop\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "      4 index 4 index moveto arc\n");
    fprintf(psfile, "      closepath\n");
    fprintf(psfile, "    gsave setrgbcolor fill grestore\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Cell fill operator:\n");
    fprintf(psfile, "%%   /xi/ /yi/ /r/ /g/ /b/ celf --> \n");
    fprintf(psfile, "/celf\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  5 3 roll \n");
    fprintf(psfile, "  %% --- r, g, b, xi, yi\n");
    fprintf(psfile, "  exch dup \n");
    fprintf(psfile, "  %% --- r, g, b, yi, xi, xi\n");
    fprintf(psfile, "  xstep mul xmin add exch 1 add xstep mul xmin add\n");
    fprintf(psfile, "  %% --- r, g, b, yi, xlo, xhi\n");
    fprintf(psfile, "  3 2 roll dup\n");
    fprintf(psfile, "  %% --- r, g, b, xlo, xhi, yi, yi\n");
    fprintf(psfile, "  ystep mul ymin add exch 1 add ystep mul ymin add\n");
    fprintf(psfile, "  %% --- r, g, b, xlo, xhi, ylo, yhi\n");
    fprintf(psfile, "  7 4 roll \n");
    fprintf(psfile, "  recf\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Triangle fill operator:\n");
    fprintf(psfile, "%%   /xa/ /ya/ /xb/ /yb/ /xc/ /yc/ /r/ /g/ /b/ trif --> \n");
    fprintf(psfile, "/trif\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    setrgbcolor\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    moveto\n");
    fprintf(psfile, "    lineto\n");
    fprintf(psfile, "    lineto\n");
    fprintf(psfile, "    closepath\n");
    fprintf(psfile, "    fill\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Segment draw operator:\n");
    fprintf(psfile, "%%   /xa/ /ya/ /xb/ /yb/ segd --> \n");
    fprintf(psfile, "/segd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "    moveto\n");
    fprintf(psfile, "    lineto\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Curve draw operator:\n");
    fprintf(psfile, "%%   /xa/ /ya/  /xb/ /yb/  /xc/ /yc/  /xd/ /yd/ arcd --> \n");
    fprintf(psfile, "/arcd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    newpath\n");
    fprintf(psfile, "      8 -2 roll moveto curveto\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Draw an X-value grid line:\n");
    fprintf(psfile, "%%   /x/ xgrd --> \n");
    fprintf(psfile, "/xgrd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "    dup ymin moveto\n");
    fprintf(psfile, "    ymax lineto\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Draw an Y-value grid line:\n");
    fprintf(psfile, "%%   /y/ ygrd --> \n");
    fprintf(psfile, "/ygrd\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "    dup xmin exch moveto\n");
    fprintf(psfile, "    xmax exch lineto\n");
    fprintf(psfile, "    stroke\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Label printing operator:\n");
    fprintf(psfile, "%%   /str/ /xa/ /ya/ /xc/ /yc/ lbsh --> \n");
    fprintf(psfile, "/lbsh\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  labelfont setfont\n");
    fprintf(psfile, "  newpath moveto\n");
    fprintf(psfile, "    %% --- str, xa, ya\n");
    fprintf(psfile, "  gsave 2 index false charpath flattenpath pathbbox grestore\n");
    fprintf(psfile, "    %% --- str, xa, ya, lox, loy, hix, hiy\n");
    fprintf(psfile, "  3 index 3 index currentpoint \n");
    fprintf(psfile, "    %% --- str, xa, ya, lox, loy, hix, hiy, lox, loy, cx, cy\n");
    fprintf(psfile, "  exch 4 1 roll exch sub\n");
    fprintf(psfile, "  3 1 roll sub exch\n");
    fprintf(psfile, "    %% --- str, xa, ya, lox, loy, hix, hiy, cx-lox, cy-loy\n");
    fprintf(psfile, "  rmoveto\n");
    fprintf(psfile, "    %% --- str, xa, ya, lox, loy, hix, hiy\n");
    fprintf(psfile, "  exch 4 1 roll exch sub \n");
    fprintf(psfile, "  3 1 roll sub exch\n");
    fprintf(psfile, "    %% --- str, xa, ya, dx, dy\n");
    fprintf(psfile, "  exch 4 1 roll mul -1 mul\n");
    fprintf(psfile, "  3 1 roll mul -1 mul exch\n");
    fprintf(psfile, "    %% --- str, -dx*xa, -dy*ya\n");
    fprintf(psfile, "  rmoveto\n");
    fprintf(psfile, "    %% --- str\n");
    fprintf(psfile, "  show\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fflush(psfile);
  }

void eps_define_caption_procs(FILE *psfile)
  {
    fprintf(psfile, "%% Operator to move to new caption line:\n");
    fprintf(psfile, "%%   nl --> \n");
    fprintf(psfile, "/nl\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  /ytext ytext dytext sub def\n");
    fprintf(psfile, "  xtext ytext moveto\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Operator to print string at CP without clipping:\n");
    fprintf(psfile, "%%   /s/ shw --> \n");
    fprintf(psfile, "/shw\n");
    fprintf(psfile, "{\n");
    fprintf(psfile, "  gsave\n");
    fprintf(psfile, "    initclip\n");
    fprintf(psfile, "    captionfont setfont show\n");
    fprintf(psfile, "  grestore\n");
    fprintf(psfile, "} def\n");
    fprintf(psfile, "\n");

    fflush(psfile);
  }

void eps_setup_page_state(FILE *psfile, int xn, int yn)
  {
    double fxn = xn;
    double fyn = yn;

    fprintf(psfile, "%% Round joints and caps:\n");
    fprintf(psfile, "1 setlinecap 1 setlinejoin\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Black thin lines:\n");
    fprintf(psfile, "0 setlinewidth 0 setgray [ ] 0 setdash\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "/xmin %f def   %% min plottable x\n", ps_hmin);
    fprintf(psfile, "/xmax %f def   %% max plottable x\n", ps_hmax);
    fprintf(psfile, "/xn %d def     %% grid cells along x axis\n", xn);
    fprintf(psfile, "/xstep %f def  %% x-size of grid cell\n", (ps_hmax-ps_hmin)/fxn);

    fprintf(psfile, "/ymin %f def   %% min plottable y\n", ps_vmin);
    fprintf(psfile, "/ymax %f def   %% max plottable y\n", ps_vmax);
    fprintf(psfile, "/yn %d def     %% grid cells along y axis\n", yn);
    fprintf(psfile, "/ystep %f def  %% y-size of grid cell\n", (ps_vmax-ps_vmin)/fyn);

    fprintf(psfile, "%% Units of measure:\n");
    fprintf(psfile, "/pt 1.0 def\n");
    fprintf(psfile, "/in pt 72.0 mul def \n");
    fprintf(psfile, "/mm pt 72.0 25.4 div mul def\n");
    fprintf(psfile, "\n");
    
    eps_set_label_font (psfile, "Courier", 8.0);    

    fprintf(psfile, "%% Set clipping path to boundary of plot area:\n");
    fprintf(psfile, "newpath\n");
    fprintf(psfile, "  xmin ymin moveto\n");
    fprintf(psfile, "  xmax ymin lineto\n");
    fprintf(psfile, "  xmax ymax lineto\n");
    fprintf(psfile, "  xmin ymax lineto\n");
    fprintf(psfile, "  xmin ymin lineto\n");
    fprintf(psfile, "clip\n");
    fprintf(psfile, "\n");

    fflush(psfile);
  }
  
void eps_setup_caption_data(FILE *psfile)
  {
    fprintf(psfile, "%% Caption text cursor:\n");
    fprintf(psfile, "/xtext xmin def\n");
    fprintf(psfile, "/ytext ymin def\n");
    fprintf(psfile, "/dytext 10 pt mul def\n");
    fprintf(psfile, "\n");

    fprintf(psfile, "%% Caption font setup:\n");
    fprintf(psfile, "/captionfont\n");
    fprintf(psfile, "  /Courier findfont\n");
    fprintf(psfile, "  dytext scalefont\n");
    fprintf(psfile, "def\n");
    fprintf(psfile, "\n");
  }
  
void eps_aux_end_page(FILE *psfile)
  {
    fprintf(psfile, "savedstate restore\n");
    fprintf(psfile, "%% Now we are back to the standard coord system.\n");
    fprintf(psfile, "\n");
    fprintf(psfile, "end %% $maindict\n");
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void eps_begin_section(FILE *psfile, char *title)
  {
    fprintf(psfile, "%%%s\n", title);
    fprintf(stderr, "[%s]\n", title);
    fflush(psfile);
  }

void eps_end_section(FILE *psfile)
  {
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void eps_add_caption(FILE *psfile, char *txt)
  {
    fprintf(psfile, "nl ");
    eps_put_text(psfile, txt, ") shw\nnl (");
    fprintf(psfile, " shw\n");
    fflush(psfile);
  }
  
void eps_set_label_font(FILE *psfile, char *font, float size)
  {
    fprintf(psfile, "%% Label font setup:\n");
    fprintf(psfile, "/labelfont\n");
    fprintf(psfile, "  /%s findfont %.3f pt mul scalefont\n", font, size);
    fprintf(psfile, "def\n");
    fprintf(psfile, "\n");
    fflush(psfile);
  }
  
void eps_put_label(
    FILE *psfile, 
    char *text, 
    double x, double y, 
    float xalign, float yalign
  )
  {
    double psx = ps_hmin + ps_xscale * (x - ps_xmin);
    double psy = ps_vmin + ps_yscale * (y - ps_ymin);
    eps_put_text(psfile, text, "\\267");
    fprintf(psfile, " ");
    fprintf(psfile, "  %5.3f %5.3f  %6.1f %6.1f  lbsh\n",
      xalign, yalign, psx, psy
    );
    fflush(psfile);
  }

void eps_draw_frame (FILE *psfile)
  {
    eps_begin_section(psfile, "Draw frame around plot area");
    fprintf(psfile, "gsave\n");
    fprintf(psfile, "%% Assumes xmax, xmin, ymax, ymin are defined.\n");
    fprintf(psfile, "  initclip\n");
    fprintf(psfile, "  newpath\n");
    fprintf(psfile, "  xmin ymin moveto\n");
    fprintf(psfile, "  xmax ymin lineto\n");
    fprintf(psfile, "  xmax ymax lineto\n");
    fprintf(psfile, "  xmin ymax lineto\n");
    fprintf(psfile, "  xmin ymin lineto\n");
    fprintf(psfile, "  closepath stroke\n");
    fprintf(psfile, "grestore\n");
    eps_end_section(psfile);
  }

void eps_set_pen(
    FILE *psfile,
    double r, double g, double b,
    double width,
    double dashlength,
    double dashspace
  )
  {
    fprintf(psfile, "%5.3f %5.3f %5.3f setrgbcolor\n", r, g, b);
    fprintf(psfile, "mm %.3f mul setlinewidth\n", width);
    if ((dashlength == 0.0) | (dashspace == 0.0))
      { fprintf(psfile, "[ ] 0 setdash\n"); }
    else
      { fprintf(psfile,
          "[ %.3f mm mul %.3f mm mul ] 0 setdash\n",
          dashlength, dashspace
        );
      }
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void eps_draw_segment(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb
  )
  {
    double psxa = ps_hmin + ps_xscale * (xa - ps_xmin);
    double psya = ps_vmin + ps_yscale * (ya - ps_ymin);
    double psxb = ps_hmin + ps_xscale * (xb - ps_xmin);
    double psyb = ps_vmin + ps_yscale * (yb - ps_ymin);
    fprintf(psfile,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxa, psya, psxb, psyb
    );
    fflush(psfile);
  }

void eps_draw_curve(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double xd, double yd
  )
  {
    double psxa = ps_hmin + ps_xscale * (xa - ps_xmin);
    double psya = ps_vmin + ps_yscale * (ya - ps_ymin);
    double psxb = ps_hmin + ps_xscale * (xb - ps_xmin);
    double psyb = ps_vmin + ps_yscale * (yb - ps_ymin);
    double psxc = ps_hmin + ps_xscale * (xc - ps_xmin);
    double psyc = ps_vmin + ps_yscale * (yc - ps_ymin);
    double psxd = ps_hmin + ps_xscale * (xd - ps_xmin);
    double psyd = ps_vmin + ps_yscale * (yd - ps_ymin);
    fprintf(psfile, "%6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f arcd\n",
      psxa, psya, psxb, psyb, psxc, psyc, psxd, psyd
    );
    fflush(psfile);
  }
  
void eps_aux_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double r, double g, double b,
    char *operator
  )
  {
    double psxlo = ps_hmin + ps_xscale * (xlo - ps_xmin);
    double psxhi = ps_hmin + ps_xscale * (xhi - ps_xmin);
    double psylo = ps_vmin + ps_yscale * (ylo - ps_ymin);
    double psyhi = ps_vmin + ps_yscale * (yhi - ps_ymin);
    fprintf(psfile, "%6.1f %6.1f  %6.1f %6.1f",
      psxlo, psxhi, psylo, psyhi
    );
    if (r >= 0.0) fprintf(psfile, "  %5.3f %5.3f %5.3f", r, g, b);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
  }

void eps_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi
  )
  {
    eps_aux_rectangle(psfile, xlo, xhi, ylo, yhi, -1.0, -1.0, -1.0, "recd");
  }

void eps_fill_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double r, double g, double b
  )
  {
    eps_aux_rectangle(psfile, xlo, xhi, ylo, yhi, r, g, b, "recf");
  }

void eps_fill_and_draw_rectangle(
    FILE *psfile,
    double xlo, double xhi,
    double ylo, double yhi,
    double r, double g, double b
  )
  {
    eps_aux_rectangle(psfile, xlo, xhi, ylo, yhi, r, g, b, "recfd");
  }

void eps_aux_polygon(
    FILE *psfile,
    double x[], double y[],
    int npoints,
    double r, double g, double b,
    char *operator
  )
  {
    int i;
    if (npoints>6) fprintf(psfile, "\n");
    for (i=0; i<npoints; i++)
      {
        double psxi = ps_hmin + ps_xscale * (x[i] - ps_xmin);
        double psyi = ps_vmin + ps_yscale * (y[i] - ps_ymin);
        if ((i % 6) == 0) 
          { if (i>0) fprintf(psfile, "\n"); }
        else
          { fprintf(psfile, "  "); }
        fprintf(psfile, "%6.1f %6.1f", psxi, psyi);
      }
    if (npoints>6) fprintf(psfile, "\n");
    fprintf(psfile, "  %d", npoints);
    if (r >= 0.0) fprintf(psfile, "  %5.3f %5.3f %5.3f", r, g, b);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
  }

void eps_draw_polygon(
    FILE *psfile,
    double x[], double y[],
    int npoints
  )
  {
    eps_aux_polygon(psfile, x, y, npoints, -1.0, -1.0, -1.0, "pold");
  }
    
void eps_fill_polygon(
    FILE *psfile,
    double x[], double y[],
    int npoints,
    double r, double g, double b
  )
  {
    eps_aux_polygon(psfile, x, y, npoints, r, g, b, "polf");
  }
    
void eps_fill_and_draw_polygon(
    FILE *psfile,
    double x[], double y[],
    int npoints,
    double r, double g, double b
  )
  {
    eps_aux_polygon(psfile, x, y, npoints, r, g, b, "polfd");
  }
    
void eps_aux_circle(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b,
    char *operator
  )
  {
    double psxc = ps_hmin + ps_xscale * (xc - ps_xmin);
    double psyc = ps_vmin + ps_yscale * (yc - ps_ymin);
    double psradius = ps_yscale * radius;
    fprintf(psfile, "%6.1f %6.1f  %6.1f", psxc, psyc, psradius); 
    if (r >= 0.0) fprintf(psfile, "  %5.3f %5.3f %5.3f", r, g, b);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
  }

void eps_fill_circle(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  )
  {
    eps_aux_circle(psfile, xc, yc, radius, r, g, b, "cirf");
  }
  
void eps_draw_circle(
    FILE *psfile,
    double xc, double yc, double radius
  )
  {
    eps_aux_circle(psfile, xc, yc, radius, -1.0, -1.0, -1.0, "cird");
  }
  
void eps_fill_and_draw_circle(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  )
  {
    eps_aux_circle(psfile, xc, yc, radius, r, g, b, "cirfd");
  }
  
void eps_aux_dot(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b,
    char *operator
  )
  {
    double psxc = ps_hmin + ps_xscale * (xc - ps_xmin);
    double psyc = ps_vmin + ps_yscale * (yc - ps_ymin);
    fprintf(psfile, "%6.1f %6.1f  %6.1f mm mul", psxc, psyc, radius); 
    if (r >= 0.0) fprintf(psfile, "  %5.3f %5.3f %5.3f", r, g, b);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
  }

void eps_fill_dot(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  )
  {
    eps_aux_dot(psfile, xc, yc, radius, r, g, b, "cirf");
  }
  
void eps_draw_dot(
    FILE *psfile,
    double xc, double yc, double radius
  )
  {
    eps_aux_dot(psfile, xc, yc, radius, -1.0, -1.0, -1.0, "cird");
  }
  
void eps_fill_and_draw_dot(
    FILE *psfile,
    double xc, double yc, double radius,
    double r, double g, double b
  )
  {
    eps_aux_dot(psfile, xc, yc, radius, r, g, b, "cirfd");
  }
  
void eps_aux_lune(
    FILE *psfile,
    double xc, double yc, double radius, double tilt,
    double r, double g, double b,
    char *operator
  )
  {
    double psxc = ps_hmin + ps_xscale * (xc - ps_xmin);
    double psyc = ps_vmin + ps_yscale * (yc - ps_ymin);
    double psradius = ps_yscale * radius;
    double pstilt = tilt * 180.0 / M_PI;
    fprintf(psfile, "%6.1f %6.1f  %6.1f %6.2f",
      psxc, psyc, psradius, pstilt); 
    if (r >= 0.0) fprintf(psfile, "  %5.3f %5.3f %5.3f", r, g, b);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
  }

void eps_fill_and_draw_lune(
    FILE *psfile,
    double xc, double yc, double radius, double tilt,
    double r, double g, double b
  )
  {
    eps_aux_lune(psfile, xc, yc, radius, tilt, r, g, b, "lunefd");
  }
  
void eps_aux_slice(
    FILE *psfile,
    double xc, double yc, double radius, double start, double stop,
    double r, double g, double b,
    char *operator
  )
  {
    double psxc = ps_hmin + ps_xscale * (xc - ps_xmin);
    double psyc = ps_vmin + ps_yscale * (yc - ps_ymin);
    double psradius = ps_yscale * radius;
    double psstart = start * 180.0 / M_PI;
    double psstop  = stop  * 180.0 / M_PI;
    fprintf(psfile, "%6.1f %6.1f  %6.1f  %6.2f %6.2f",
      psxc, psyc, psradius, psstart, psstop); 
    if (r >= 0.0) fprintf(psfile, "  %5.3f %5.3f %5.3f", r, g, b);
    fprintf(psfile, " %s\n", operator);
    fflush(psfile);
  }

void eps_fill_and_draw_slice(
    FILE *psfile,
    double xc, double yc, double radius, double start, double stop,
    double r, double g, double b
  )
  {
    eps_aux_slice(psfile, xc, yc, radius, start, stop, r, g, b, "slicefd");
  }
  
void eps_fill_triangle(
    FILE *psfile,
    double xa, double ya,
    double xb, double yb,
    double xc, double yc,
    double r, double g, double b
  )
  {
    double psxa = ps_hmin + ps_xscale * (xa - ps_xmin);
    double psya = ps_vmin + ps_yscale * (ya - ps_ymin);
    double psxb = ps_hmin + ps_xscale * (xb - ps_xmin);
    double psyb = ps_vmin + ps_yscale * (yb - ps_ymin);
    double psxc = ps_hmin + ps_xscale * (xc - ps_xmin);
    double psyc = ps_vmin + ps_yscale * (yc - ps_ymin);
    fprintf(psfile,
      "%6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f  %5.3f %5.3f %5.3f trif\n",
            psxa, psya, psxb, psyb, psxc, psyc, r, g, b
    );
    fflush(psfile);
  }

void eps_fill_grid_cell(FILE *psfile, int xi, int yi, double r, double g, double b)
  {
    fprintf(psfile, "%3d %3d  %5.3f %5.3f %5.3f celf\n", xi, yi, r, g, b);
    fflush(psfile);
  }

void eps_draw_coord_line (FILE *psfile, char axis, double coord)
  {
    double pscoord;
    char psaxis;
    if ((axis == 'x') || (axis == 'X'))
      { pscoord = ps_hmin + ps_xscale * (coord - ps_xmin); psaxis = 'x'; }
    else if ((axis == 'y') || (axis == 'Y'))
      { pscoord = ps_vmin + ps_yscale * (coord - ps_ymin); psaxis = 'y'; }
    else
      { error("eps_draw_coord-line: invalid axis"); }
    fprintf(psfile, "%6.1f %cgrd\n", pscoord, psaxis);
  }

void eps_draw_grid_lines(FILE *psfile)
  {
    fprintf(psfile, "%% Grid lines:\n");
    fprintf(psfile, "gsave\n");
    fprintf(psfile, "  initclip\n");
    fprintf(psfile, "  0 1 xn {\n");
    fprintf(psfile, "    xstep mul xmin add xgrd\n");
    fprintf(psfile, "  } for\n");
    fprintf(psfile, "  0 1 yn {\n");
    fprintf(psfile, "    ystep mul ymin add ygrd\n");
    fprintf(psfile, "  } for\n");
    fprintf(psfile, "grestore\n");
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void eps_put_text(FILE *psfile, char *text, char *newline)
  {
    char *p;
    putc('(', psfile);
    for (p = text; *p != 0; p++)
      {
        if (*p == '\n')
          { fprintf(psfile, "%s", newline); }
        else if (*p == '(')
          { putc('\\', psfile); putc('(', psfile); }
        else if (*p == ')')
          { putc('\\', psfile); putc(')', psfile); }
        else if (*p == '\t')
          { putc(' ', psfile); putc(' ', psfile); }
        else if (*p == '\\')
          { putc('\\', psfile); putc('\\', psfile); }
        else if ((*p < ' ') || (*p > '~'))
          { fprintf(psfile, "\\%03o", *p); }
        else
          { putc(*p, psfile); }
      }
    fprintf(psfile, ")");
  }
