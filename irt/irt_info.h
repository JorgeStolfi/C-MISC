#ifndef irt_info_H
#define irt_info_H

/* Last edited on 2023-02-22 11:48:18 by stolfi */
/* Documentation of the Interval Ray Tracer. */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  {DIR} {SCENE_NAME}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Produces an image of a 3D scene by tay-tracing based on\n" \
  " Interval Arithmetic (IA) or Affine Arithmetic (AA).\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  ???.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -foo ??? \n" \
  "    This optional argument specifies ???\n" \
  "\n" \
  "  -bar ??? \n" \
  "    This mandatory argument specifies ???\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  ??? (1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  A first version of this program was written in 1993 by Marcus V. A. Andrade, " \
  " vaguely inspired on Roman Kuchkuda's KRT ray tracer.  It was substantially" \
  " modified by J. Stolfi in 1994--1995 and again in 2007 and again in 2008.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " irt_main_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  

#endif
