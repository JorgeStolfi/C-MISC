#define PROG_NAME "ioplot"
#define PROG_DESC "create a sliced-ham in-out plot in Encapsulated Postscript"
#define PROG_VERS "2013-10-27"
/* Last edited on 2024-12-21 11:56:05 by stolfi */

#define PROG_COPYRIGHT "© 2005  State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -figSize {MMWIDTH} {MMHEIGHT} ] \\\n" \
  "  [ -plotWidth {MM} ] \\\n" \
  "  [ -scale {VALPERMM} ] \\\n" \
  "  [ -leafGap {MM} ] \\\n" \
  "  [ -labelGap {MM} ] \\\n" \
  "  [ -showNodes ] \\\n" \
  "  [ -font {NAME} {PTSIZE} ] [ -fmt {VALFMT} ] \\\n" \
  "  [ -fillColor {RF} {GF} {BF} ] \\\n" \
  "  [ -drawColor {RF} {GF} {BF} ] \\\n" \
  "  [ -textColor {RF} {GF} {BF} ] \\\n" \
  "  [ -depth {POSDEPTH} {NEGDEPTH} ] \\\n" \
  "  [ -simplify ] \\\n" \
  "  [ -debug ] \\\n" \
  "  < {INFILE}.txt  > {OUTFILE}.eps"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads a file with positive and negative entries, e.g. the" \
  " items of a budget or financial report.  Outputs a" \
  " multilevel \"sliced ham\" diagram showing positive" \
  " entries on the left side, negative entries on the" \
  " right side, each with its corresponding label.\n" \
  "\n" \
  "  Each entry of the diagram is defined by one line of the input" \
  " file, in the format \"{INDENT} {VALUE} {LABEL}\".  Here the" \
  " {INDENT} field consists of zero or more '>' characters, {VALUE} is either '*' or" \
  " a decimal number (possibly signed and/or frctional), and {LABEL} is any text. There should be" \
  " at least one space between {VALUE} and {LABEL}; additional spaces may also be" \
  " inserted before, between, or after the other fields.\n" \
  "\n" \
  "  Entries with {VALUE} '+0' or '-0' are treated as positive" \
  " or negative, respectively; unsigned" \
  " nonzero {VALUE}s are assumed positive, and entries with unsigned" \
  " zero value are ignored.\n" \
  "\n" \
  "  The optional {INDENT} fields organize the entries into" \
  " two trees, one with all the positive items, and one with" \
  " all the negative items.  The depth of each entry in the plotted" \
  " tree is 1 plus the number of leading '>'s.  An entry depth {N} is" \
  " assumed to be a sub-item of the most recent entry with depth {N-1}.\n" \
  "\n" \
  "  The leaves of each tree must have a numeric {VALUE} field, while" \
  " the internal nodes must have {VALUE} == '*', standing for the sum of all" \
  " leaves in its subtree.  Each internal node will be assigned" \
  " to the positive tree, to the negative tree, or to both" \
  " trees, depending on the signs of the {VALUE}s of its" \
  " leaves.\n" \
  "\n" \
  "  In the output diagram, the entries in each tree are plotted" \
  " bottom to top.  Each entry is plotted as a rectangle whose" \
  " width is fixed, and whose height is" \
  " proportional to its absolute {VALUE}.  Nodes at the same" \
  " depth are plotted on" \
  " the same column of the diagram.  The root of each tree is" \
  " a single rectangle near the middle of the diagram.  Each" \
  " node is connected to its parent by a band as tall as itself.\n" \
  "\n" \
  "OPTIONS\n" \
  "  Except for the font size, all dimensions are in mm.\n" \
  "\n" \
  "  -figSize {MMWIDTH} {MMHEIGHT}\n" \
  "    Specifies the size of the whole figure," \
  "including labels and margins.\n" \
  "\n" \
  "  -plotWidth {MM}\n" \
  "    Specifies the total width of the sliced-ham diagram" \
  " proper, excluding the labels.  If this parameter is not" \
  " given, the plot uses the central 1/2 of the total figure width," \
  " minus margins.\n" \
  "\n" \
  "  -scale {VALPERMM}\n" \
  "    Specifies the vertical scale of the plot, namely" \
  " what {VALUE} is to be plotted as 1 mm.  If this parameter" \
  " is not given, the vertical scale is set so that the plot fits in the total" \
  " figure height, minus the margins.\n" \
  "\n" \
  "  -leafGap {MM}\n" \
  "    Specifies the vertical space between the free ends of the slices.\n" \
  "\n" \
  "  -labelGap {MM}\n" \
  "    Defines the horizontal spacing between the plot and the labels.\n" \
  "\n" \
  "  -showNodes\n" \
  "    If specified, highlights each leaf or internal node by drawing" \
  " vertical lines on both sides of it.  Otherwise draws the outline" \
  " only of the whole tree, so that each node is fused with the color bands" \
  " that connect it to its parent and children.\n" \
  "\n" \
  "  -fmt {VALFMT}\n" \
  "    Specifies the format string (as in {printf(3C)}) that the program" \
  " will use to print the data values alongside the plot.\n" \
  "\n" \
  "  -font {NAME} {PTSIZE}\n" \
  "    Specifies the name and size (in pt) of" \
  " the font that will be used for values and labels.\n" \
  "\n" \
  "  -fillColor {RF} {GF} {BF}\n" \
  "  -drawColor {RF} {GF} {BF}\n" \
  "  -textColor {RF} {GF} {BF}\n" \
  "    Specify the colors for the tree interior, tree outline," \
  " and labels, respectively. The arguments are the R,G,B " \
  " intensities, as fractions between 0 0 0 (black) and" \
  " 1 1 1 (white).  If any component is negative, " \
  " the corresponding layer is omitted from the figure.\n" \
  "\n" \
  "  -depth {POSDEPTH} {NEGDEPTH}\n" \
  "    Specifies the depth at which each tree should be truncated.  Any entries" \
  " below that depth are omitted. Any entry at that depth is plotted as a" \
  " leaf node, with value equal to the sum of all the omitted leaves" \
  " under it.  In particular, depth 0 results in a single slice" \
  " showing the total value; depth -1 omits the tree altogether; and" \
  " depth \"all\" shows whe whole tree.\n" \
  "\n" \
  "  -simplify\n" \
  "    Specifies that internal nodes with only one child should be omitted.\n" \
  "\n" \
  "  -debug\n" \
  "    Requests a printout (to {stderr}) showing the input entries and also" \
  " the computed values of internal nodes.\n" \
  "\n" \
  "EXAMPLE\n" \
  "  Here is a sample input:\n" \
  "\n" \
  "    +250 Initial Cash Bal\n" \
  "    -35 Final Cash Bal\n" \
  "    +50 Sales\n" \
  "    +30 Gifts Received\n" \
  "    -80 Rent\n" \
  "    * Personel\n" \
  "    > -20 Lawyers\n" \
  "    > -10 Janitors\n" \
  "    * Services\n" \
  "    > +5 Consulting\n" \
  "    > -30 Royalties\n" \
  "    * Plots\n" \
  "    > * Color plots\n" \
  "    > > -20 Blue plots\n" \
  "    > > -30 Green plots\n" \
  "    > * B&W plots\n" \
  "    > > -20 Piecharts\n" \
  "    > > -20 Bargraphs\n" \
  "    > -30 Murky plots\n" \
  "    > +15 Plot refunds\n" \
  "    > +25 Plot sales\n" \
  "    -70 Food\n" \
  "    -10 Toner\n" \
  "\n" \
  "  Here are the two trees that will be plotted from this input:\n" \
  "\n" \
  "    Positive tree                      Negative tree\n" \
  "    --------------------------------   ---------------------------------\n" \
  "    Total . . . . . . .  375           Total . . . . . . 375            \n" \
  "    + Initial Cash Bal  . .  250       + Final Cash Bal  . .  35        \n" \
  "    + Sales . . . . . . . . . 50       + Rent  . . . . . . .  80        \n" \
  "    + Gifts Received  . . . . 30       + Personel  . . . . .  30        \n" \
  "    |                                  | + Lawyers . . . . . . .  20    \n" \
  "    + Services  . . . . . . .  5       | + Janitors  . . . . . .  10    \n" \
  "    | + Consulting  . . . . . . .  5   |                                \n" \
  "    |                                  + Services  . . . . .  30        \n" \
  "    + Plots . . . . . . . .   40       | + Royalties . . . . . .  30    \n" \
  "    . + Plot refunds  . . . . . . 15   |                                \n" \
  "    . + Plot sales  . . . . . . . 25   + Plots . . . . . .   130        \n" \
  "                                       | + Color plots . . . . .  50    \n" \
  "                                       | | + Blue plots. . . . . . .  20\n" \
  "                                       | | + Green plots . . . . . .  30\n" \
  "                                       | |                              \n" \
  "                                       | + B&W plots . . . . . .  40    \n" \
  "                                       | | + Piecharts . . . . . . .  20\n" \
  "                                       | | + Bargraphs . . . . . . .  20\n" \
  "                                       | |                              \n" \
  "                                       | + Murky plots . . . . .  30    \n" \
  "                                       |                                \n" \
  "                                       + Food  . . . . . . .  70        \n" \
  "                                       + Toner . . . . . . .  10        \n" \
  "\n" \
  "  Note that the final cash balance, which is actually +35, is here entered" \
  " as -35 so that it appears on the negative tree (and so that the totals" \
  " would match).\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  gnuplot(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in March 2005 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2009-10-31:\n" \
  "    Added \"--info\" option. Added \"-showNodes\" option.\n" \
  " Renamed \"-plotColor\", \"-lineColor\", \"fontColor\"" \
  " to \"-fillColor\", \"-drawColor\", \"textColor\". [J.Stolfi].\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <math.h>

#include <epswr.h>
#include <bool.h>
#include <sign.h>
#include <frgb.h>
#include <affirm.h>
#include <jsstring.h>
#include <argparser.h>

typedef struct options_t
  { double figSizeX;   /* Width of figure (mm). */
    double figSizeY;   /* Height of figure (mm). */
    double scale;      /* Vertical scale (valunits per mm); 0 to auto-scale. */
    double plotWidth;  /* Total X-width of sliced ham (mm). */
    double leafGap;    /* Y-gap between slices (mm). */
    double labelGap;   /* X-space between ham and labels (mm). */
    bool_t showNodes;  /* TRUE to draw both vertical sides of every node. */
    char *fmt;         /* Format for values, or NULL. */
    char *fontName;    /* Label font name, e.g. "TimesRoman", "Helvetica", etc.. */
    double fontSize;   /* Label font size (pt). */
    frgb_t textColor;  /* Color for value and label text. */
    frgb_t fillColor;  /* Color for tree plot interior. */
    frgb_t drawColor;  /* Color for tree plot outline. */
    int32_t depth_pos;     /* Max plot depth for positive tree; -2 for input depth. */
    int32_t depth_neg;     /* Max plot depth for negative tree; -2 for input depth. */
    bool_t simplify;   /* TRUE deletes internal nodes with only one child. */
    bool_t debug;      /* TRUE prints entries with subtree totals. */
  } options_t;

#define DEFAULT_figSizeX (150.0)
#define DEFAULT_figSizeY (150.0)
#define DEFAULT_leafGap (2.0)
#define DEFAULT_fillColor ((frgb_t){{ 1.000f, 0.800f, 0.700f }})
#define DEFAULT_drawColor ((frgb_t){{ 0.200f, 0.000f, 0.000f }})
#define DEFAULT_fontName "Courier"
#define DEFAULT_fontSize (8.0)
#define DEFAULT_textColor ((frgb_t){{ 0.000f, 0.000f, 0.000f }})
#define DEFAULT_fmt "%.2f"
#define DEFAULT_labelGap (2.0)
#define DEFAULT_showNodes (TRUE)
  /* Basic default values for command line options. Other
    defaults are derived from these or from given values. */

#define DEFAULT_plotWidthFraction (1.0/3.0)
  /* Default for {plotWidth} relative to {figSizeX}. */

/*
  DATA TREES */

typedef struct item_t
  { double val;  /* Sum of absolute value of all leaf items under this node. */
    char *name;  /* Label to print for entry. */
    struct item_t *sub;   /* First child, or NULL if leaf. */
    struct item_t *next;  /* Next sibling. */
  } item_t;
  /*
    Input data lines of same sign are represented as nodes of a tree;
    each {item_t} record represents one data item. If there are any
    nonzero items in the tree, the root of the tree is a dummy item
    with no siblings and name "TOTAL", whose children are the
    top-level items. */

#define MAXDEPTH 255
  /* Maximum allowed depth (absurdly high). */

bool_t read_item(FILE *rd, int32_t *lineP, /*OUT*/ item_t **itP, int32_t *dpnP, sign_t *sideP, bool_t *leafP );
  /* Tries to read an entry from {stdin}, skipping blank and comment
    lines. If it succeeds, stores the address of the entry's record in
    {*itP}, its depth (number of leading '>'s) in {*dpnP}, its sign in
    {*sideP}, sets {*leafP} to TRUE iff the node is a leaf, and returns
    {TRUE}. Returns {FALSE} on EOF. In any case, increments the
    {*line} counter appropriately. */

item_t *make_root_item( void );
  /* Creates an item record with zero value, no children
  or siblings, and a fresh copy of "TOTAL" for {name}. */

item_t *copy_item ( item_t *it );
  /* Makes a copy of item {it}, including its {name} string. */

void read_data ( FILE *rd, /*OUT*/ item_t **posP, item_t **negP );
  /* Reads entries from a data file, and builds separate trees for the
    positive and negative ones. Omits leaves with unsigned zero
    {VALUE}, and internal nodes with no children. Returns in {*posP}; and analogously
    for {*negP}. */

void free_tree ( item_t *it );
  /* Reclaims all nodes of the tree rooted at {it},
    including the {name} strings. */

item_t *trim_tree ( item_t *it, int32_t dp );
  /* Removes all levels of the tree {it} beyond level {dp}.
    In particular, if {dp = 0} leaves only the root node.
    Any internal nodes at level {dp} becomes a leaf. */

item_t *simplify_tree ( item_t *it );
  /* Removes any non-leaf nodes that have only one child. */

void print_tree ( FILE *wr, char *fmt, item_t *it );
  /* Prints the tree rooted at {it} to {wr}. */

item_t *reverse_siblings( item_t *it );
  /* Destructively reverse the order of {it} and its siblings
    Does not affect the proper descendants of those nodes. */

/*
  TREE STATISTICS */

void summarize_tree ( item_t *it, int32_t *nvP, int32_t *dpmaxP, double *totP );
  /* Computes the number of leaves {*nvP}, the maximum leaf depth {*dpmaxP},
    and the sum of all leaf values {*totP}. */

/*
  PLOTTING */

#define REL_NODE_WD (0.50)
  /* Width of node rectangle relative to level width. */

double pick_scale ( int32_t nv, double tot, double gap, double ysize );
  /* Selects the scale for the vertical axis, given the absolute total
    {tot} and the number of leaf items {nv} on one of the trees.
    Assumes that the figure is {ysize} mm tall (excluding margins) and
    that the slices are spread out {gap} mm apart at their free end. */

void pick_abscissas( options_t *o, int32_t dpos, int32_t dneg, double *xmidP, double *xposP, double *xnegP);
  /* Selects the abscissa {*xmidP} for the root end of both plots, and
    the abscissas {*xposP,*xnegP} of the leaf tips for the positive tree
    (on the left side) and negative tree (on the right side). Assumes that
    their depths are {dneg} and {dpos}, respectively. */


epswr_figure_t *init_plot_file
  ( char *dir, 
    char *prefix, 
    char *name, 
    int32_t seq, 
    char *suffix, 
    options_t *o
  );
  /* Create an Encapsulated Postscript plotting stream,
    that writes to a file called "{dir}/{prefix}_{name}_{NNNNN}_{suffix}.eps".
    Any name that are {NULL} or {-1} are omitted, with the relevant separators.
    The scale {o->scale} must be properly set. */

void fini_plot_file ( epswr_figure_t *eps, options_t *o );
  /* Terminates the Postscript stream. */

void plot_slices
  ( epswr_figure_t *eps,  /* Postscript figure stream. */
    options_t *o,  /* Command line options. */
    double xbas,   /* X coordinate of uncut edge of ham. */
    double xtip,   /* X coordinate of cut edge of ham. */
    item_t *it,    /* The root node, or NULL. */
    int32_t dpmax,     /* Treat any nodes at this depth as leaves. */
    sign_t side    /* Side of plot. */
  );
  /* Plots one half of the sliced-ham diagram: {side = +1} means the
    left half, {side = -1} the right half. Assumes that the total is
    {tot} and the maximum leaf depth is {dpmax}. The item list
    at each level is plotted from bottom to top. */

void plot_branch_polygonal
  ( epswr_figure_t *eps,  /* Postscript figure stream. */
    options_t *o,  /* Command line options. */
    double xjoin,  /* X coordinate of joined (parent) side of branch. */
    double xopen,  /* X coordinate of open (child) side of branch. */
    double yjoin,  /* Bottom Y of branch, joined (parent) side. */
    double yopen,  /* Bottom Y of branch, open (child) side. */
    double dy,     /* Y width of branch. */
    int32_t pass       /* Plotting pass. */
  );
  /* Plots the band connecting an internal node to a child node.
    The X range spans from {xjoin} (parent side) to {xopen} (child side).
    Note that these two abscissas may be in any order.
    The Y range is {[yjoin _ yjoin+dy]} at the join (parent) sie,
    and {[yopen _ yopen+dy]} at the open (child) side. */

void plot_branch_curved
  ( epswr_figure_t *eps,  /* Postscript figure stream. */
    options_t *o,  /* Command line options. */
    double xjoin,  /* X coordinate of joined (parent) side of branch. */
    double xopen,  /* X coordinate of open (child) side of branch. */
    double yjoin,  /* Bottom Y of branch, joined (parent) side. */
    double yopen,  /* Bottom Y of branch, open (child) side. */
    double dy,     /* Y width of branch. */
    int32_t pass       /* Plotting pass. */
  );
  /* Same as plot_band_polygonal, but uses Bézier arcs intead of broken lines. */

/* OTHER PROTOTYPES */

int32_t main ( int32_t argc, char **argv );

options_t *get_options ( int32_t argc, char **argv );
frgb_t parse_frgb ( argparser_t *pp );
int32_t parse_depth ( argparser_t *pp );
void data_error ( int32_t line, char *msg );

/* IMPLEMENTATIONS */

#define verbose FALSE

int32_t main (int32_t argc, char **argv)
  {
    /* Parse command line arguments: */
    options_t *o = get_options(argc, argv);

    item_t *pos; item_t *neg;
    read_data(stdin, &pos, &neg);

    /* Truncate the two trees: */
    if (o->depth_pos != MAXDEPTH) { pos = trim_tree(pos, o->depth_pos); }
    if (o->depth_neg != MAXDEPTH) { neg = trim_tree(neg, o->depth_neg); }

    if (o->simplify)
      { /* Remove trivial internal nodes: */
        pos = simplify_tree(pos);
        neg = simplify_tree(neg);
      }

    if (o->debug)
      { fprintf(stderr, "positive tree:\n"); print_tree(stderr, o->fmt, pos);
        fprintf(stderr, "\n");
        fprintf(stderr, "negative tree:\n"); print_tree(stderr, o->fmt, neg);
      }

    /* Summarize trees: */
    int32_t dpos, dneg;     /* Maximum depth of any node, or -1 if empty. */
    double tpos, tneg;  /* Sum of leaf values. */
    int32_t npos, nneg;     /* Number of leaves in each tree. */
    summarize_tree(pos, &npos, &dpos, &tpos);
    summarize_tree(neg, &nneg, &dneg, &tneg);

    /* If the user specified a depth, force that depth: */
    if (o->depth_pos != MAXDEPTH) { dpos = o->depth_pos; }
    if (o->depth_neg != MAXDEPTH) { dneg = o->depth_neg; }

    /* Compute scale if not specified: */
    if (o->scale <= 0)
      { double pos_scale = pick_scale(npos, tpos, o->leafGap, o->figSizeY);
        double neg_scale = pick_scale(nneg, tneg, o->leafGap, o->figSizeY);
        o->scale = (pos_scale > neg_scale ? pos_scale : neg_scale);
      }
    if (o->debug)
      { /* Print scale: */
        fprintf(stderr, "vertical scale = %.6f/mm\n", o->scale);
      }

    /* Compute the root abscissa {xmid} and the leaf abscissas {xneg,xpos}: */
    double xmid, xpos, xneg;
    pick_abscissas(o, dpos, dneg, &xmid, &xpos, &xneg);

    epswr_figure_t *eps = init_plot_file(NULL, NULL, NULL, -1, NULL, o);
    plot_slices(eps, o, xmid, xpos, pos, dpos, +1);
    plot_slices(eps, o, xmid, xneg, neg, dneg, -1);
    fini_plot_file(eps, o);
    return(0);
  }

void read_data ( FILE *rd, /*OUT*/ item_t **posP, item_t **negP )
  {
    /* Positive and negative stacks: */
    int32_t pos_dps;
    item_t *pos_stk[MAXDEPTH+1];
    double pos_tot[MAXDEPTH+1];

    int32_t neg_dps;
    item_t *neg_stk[MAXDEPTH+1];
    double neg_tot[MAXDEPTH+1];

    /* During the procedure, {pos_dps} is the depth of the last
      positive item read,  and {pos_stk[k]} points to
      the last positive item of depth {k}, for each {k} in
      {0..pos_dps}. Let {it = pos_stk[k]}; the field {it->next} points
      to the /previous/ sibling of {it} (or NULL). The variable
      {pos_tot[k]} holds the the sum of all values in the leaves
      hanging from {it} and its siblings.

      Ditto for {neg_dps}, {neg_stk}, etc.. */

    int32_t line = 0; /* Current input line number. */

    auto void finish_nephews ( int32_t *dpsP, item_t *stk[], double tot[], int32_t dpn );

    auto void add_item ( int32_t *dpsP, item_t *stk[], double tot[], item_t *it, int32_t dpn, bool_t leaf );
      /* Adds to the stack {stk[0..*dps],tot[0..*dps]} an item {it} with depth {dpn},
        which may be leaf or not according to the {leaf} parameter. */

    void finish_nephews ( int32_t *dpsP, item_t *stk[], double tot[], int32_t dpn )
      {
        int32_t dps = (*dpsP);

        /* Close all currently open subtrees deeper than {dpn}: */
        while (dps > dpn)
          {
            /* Terminates the current sibling list {stk[*dps]}. */
            item_t *its = stk[dps];
            /* Get parent node: */
            int32_t dpp = dps-1; item_t *itp = stk[dpp];
            if (its == NULL)
              { /* No children, omit parent node: */
                stk[dpp] = itp->next;
                free(itp);
              }
            else
              { /* Reverse sibling list and hang from parent: */
                itp->sub = reverse_siblings(its);
                /* Set the parent's value to be the subtree sum: */
                itp->val = tot[dps];
                tot[dpp] += itp->val;
              }
            dps = dpp;
          }

        (*dpsP) = dps;
      }

    void add_item ( int32_t *dpsP, item_t *stk[], double tot[], item_t *it, int32_t dpn, bool_t leaf )
      {
        assert(dpn >= 0);
        int32_t dps = (*dpsP);

        /* Check for proper nesting: */
        if (dpn > dps)
          { data_error(line, "parent \'*\' entry is missing"); }

        /* Close trees at higher levels, if any: */
        finish_nephews(&dps, stk, tot, dpn);
        assert(dpn == dps);

        /* Add item {it} to the current sibling list: */
        it->next = stk[dps];
        stk[dps] = it;
        it->sub = NULL;
        if (leaf)
          { /* Add leaf value to list: */
            tot[dps] += it->val;
          }
        else
          { /* Start children list: */
            int32_t dpc = dps + 1;
            assert(dpc <= MAXDEPTH);
            stk[dpc] = NULL; tot[dpc] = 0;
            /* Next item should be of depth {dpc} or less: */
            dps = dpc;
          }

        (*dpsP) = dps;
      }

    item_t *it;
    int32_t dpn;
    bool_t leaf;

    /* Initialize the stacks to accept the dummy root item: */
    pos_dps = 0; pos_stk[0] = NULL; pos_tot[0] = 0.0;
    neg_dps = 0; neg_stk[0] = NULL; neg_tot[0] = 0.0;

    /* Add the dummy root items: */
    add_item(&pos_dps, pos_stk, pos_tot, make_root_item(), 0, FALSE);
    add_item(&neg_dps, neg_stk, neg_tot, make_root_item(), 0, FALSE);

    /* Add the real items (with depth >= 1): */
    sign_t side;
    while(read_item(stdin, &line, &it, &dpn, &side, &leaf))
      { /* Add item to proper tree: */
        if (side > 0)
          { add_item(&pos_dps, pos_stk, pos_tot, it, dpn, leaf); }
        else if (side < 0)
          { add_item(&neg_dps, neg_stk, neg_tot, it, dpn, leaf); }
        else if (! leaf)
          { /* Internal node, add to both lists: */
            item_t *it_pos = it;
            item_t *it_neg = copy_item(it);
            add_item(&pos_dps, pos_stk, pos_tot, it_pos, dpn, leaf);
            add_item(&neg_dps, neg_stk, neg_tot, it_neg, dpn, leaf);
          }
        else
          { /* Ignore zero-sign leaf items. */
            free_tree(it);
          }
      }

    /* Close all positive levels, so that root is OK: */
    finish_nephews(&pos_dps, pos_stk, pos_tot, 0); (*posP) = pos_stk[0];
    finish_nephews(&neg_dps, neg_stk, neg_tot, 0); (*negP) = neg_stk[0];
  }

item_t *make_root_item(void)
  { item_t *r = malloc(sizeof(item_t));
    r->val = 0.0;
    r->name = txtcat("TOTAL", ""); /* New copy of "TOTAL" */
    r->next = NULL;
    r->sub = NULL;
    return r;
  }

item_t *copy_item ( item_t *it )
  { item_t *r = malloc(sizeof(item_t));
    r->val = it->val;
    r->name = txtcat(it->name, ""); /* New copy of {name} string. */
    r->next = it->next;
    r->sub = it->sub;
    return r;
  }

bool_t read_item(FILE *rd, int32_t *lineP, /*OUT*/ item_t **itP, int32_t *dpnP, sign_t *sideP, bool_t *leafP )
  {
    int32_t c = getc(rd);
    int32_t sgn;    /* Sign of entry. */
    double val; /* Value of entry (including sign): */
    int32_t dpn;    /* Depth of entry (starting at 1). */
    while (1)
      {
        /* Check for EOF at beginning-of-line: */
        if (c == EOF) { return FALSE; }

        /* Got one more line: */
        (*lineP)++;

        /* Skip leading blanks and indentators, counting the latter:  */
        dpn = 1;
        while (1)
          {
            /* Skip leading blanks: */
            while (
                (c == ' ') ||
                (c == '\011') ||
                (c == '\014') ||
                (c == '\015') ||
                (c == '\240')
              )
              { c = getc(rd); }

            /* If line starts with '#', skip rest of line: */
            if (c == '#')
              { do { c = getc(rd); } while ((c != '\012') && (c != EOF)); }

            if (c == '>')
              { dpn++; c = getc(rd); }
            else
              { break; }
          }

        /* What have we got now: */
        if (c == EOF)
          { data_error(*lineP, "unexpected end of file"); }
        else if (c == '\n')
          { /* Line was blank or comment, ignore: */
            if (dpn > 1)
              { data_error(*lineP, "missing value field"); }
            /* Start next line: */
            c = getc(rd);
            continue;
          }
        else if (c == '*')
          { /* Internal node (implies one extra level): */
            if (dpn > MAXDEPTH-1)
              { data_error(*lineP, "too many nesting levels"); }
            (*leafP) = FALSE;
            sgn = 0;
            val = 0.0;
            break;
          }
        else if ((c == '+') || (c == '-') || ((c >= '0') && (c <= '9')) || (c == '.'))
          { /* Leaf node: */
            (*leafP) = TRUE;
            if (dpn > MAXDEPTH)
              { data_error(*lineP, "too many nesting levels"); }
            if (c == '+')
              { sgn = +1; }
            else if (c == '-')
              { sgn = -1; }
            else
              { sgn = 0 /* (for now) */; break; }
            ungetc(c, rd);

            /* Read value from line: */
            if (fscanf(rd, "%lf", &val) != 1)
              { data_error(*lineP, "invalid number"); }

            /* Adjust sign if needed: */
            if (sgn > 0)
              { /* Explicit '+' sign: */ assert(val >= 0); }
            else if (sgn < 0)
              { /* Explicit '-' sign: */ assert(val <= 0); }
            else
              { /* No explicit sign: */ assert(val >= 0); if (val > 0) { sgn = +1; } }
            break;
          }
        else
          { data_error(*lineP, "item value invalid or missing"); }
      }

    /* Require at least one blank space: */
    c = getc(rd);
    if (c == EOF)
      { data_error(*lineP, "unexpected end of file"); }

    /* Read label from line: */
    size_t len;
    ssize_t nread = getline(&name, &len, rd);
    if (nread < 0)
      { data_error(*lineP, "unexpected end of file"); }
    if (nread == 0)
      { data_error(*lineP, "missing label"); }
    name[nread-1] = '\0';

    /* Allocate record: */
    item_t *it = (item_t *)notnull(malloc(sizeof(item_t)), "no mem");

    /* Fill and return: */
    it->val = sgn*val; /* Absolute value of {val}. */
    it->name = name;
    it->next = NULL;
    (*sideP) = sgn;
    (*itP) = it;
    (*dpnP) = dpn;
    return TRUE;
  }

void free_tree ( item_t *it )
  { while (it != NULL)
      { item_t *nxt = it->next;
        if (it->sub != NULL)
          { free_tree(it->sub); }
        if (it->name != NULL)
          { free(it->name); }
        free(it);
        it = nxt;
      }
  }

item_t *trim_tree ( item_t *it, int32_t dp )
  {
    if (it == NULL)
      { return NULL; }
    else if (dp < 0)
      { free_tree(it); return NULL; }
    else
      { item_t *t = it;
        while (t != NULL)
          { if (t->sub != NULL)
              { t->sub = trim_tree(t->sub, dp-1); }
            t = t->next;
          }
        return it;
      }
  }

item_t *simplify_tree ( item_t *it )
  {
    if (it == NULL)
      { return NULL; }
    else
      { item_t *r = NULL;
        while (it != NULL)
          { item_t *n = it->next;
            item_t *s = simplify_tree(it->sub);
            if ((s != NULL) && (s->next == NULL))
              { /* {s} is the only child of {it}: */ free(it); it = s; }
            else
              { /* {it} is leaf, or {s} has siblings: */ it->sub = s; }
            it->next = r; r = it; it = n;
          }
        return reverse_siblings(r);
      }
  }

void summarize_tree ( item_t *it, int32_t *nvP, int32_t *dpmaxP, double *totP )
  {
    int32_t nv = 0;
    int32_t dpmax = -1;
    double tot = 0.0;

    auto void do_sum ( item_t *t, int32_t dp );
    /* Summarize the subtree rooted at {t}. If not NULL,
       its items are assumed to be at level {dp}. */

    void do_sum( item_t *t, int32_t dp )
      { while (t != NULL)
          { if (t->sub == NULL)
              { /* Leaf: */
                nv++;
                tot += t->val;
                if (dp > dpmax) { dpmax = dp; }
              }
            else
              { /* Non-leaf: */
                do_sum(t->sub, dp+1);
              }
            t = t->next;
          }
      }

    do_sum(it, 0);
    (*nvP) = nv;
    (*dpmaxP) = dpmax;
    (*totP) = tot;
  }

void print_tree ( FILE *wr, char *fmt, item_t *it )
  {

    auto void do_print ( item_t *t, int32_t dp );
    /* Print the subtree rooted at {t}. If not NULL,
       its items are assumed to be at level {dp}. */

    void do_print( item_t *t, int32_t dp )
      { while (t != NULL)
          { for (uint32_t k = 0;  k < dp; k++) { fprintf(wr, "> "); }
            if (t->sub == NULL)
              { fprintf(wr, fmt, t->val); 
                fprintf(wr, " %s\n", t->name);
              }
            else
              { fprintf(wr, "*(=");
                fprintf(wr, fmt, t->val); 
                fprintf(wr, ") %s\n", t->name);
                do_print(t->sub, dp+1);
              }
            t = t->next;
          }
      }

    do_print(it, 0);
  }

double pick_scale ( int32_t nv, double tot, double gap, double ysize )
  {
    /* Compute vertical extent of figure, minus all gaps: */
    double yfree = ysize - nv*gap;
    /* Fit slices if possible, else fit unsliced ham in half the height: */
    double scale = tot / (yfree > 0 ? yfree : ysize/2);
    return scale;
  }

void pick_abscissas( options_t *o, int32_t dpos, int32_t dneg, double *xmidP, double *xposP, double *xnegP)
  {

    /* Center the plot within the figure: */
    double hamwd = o->plotWidth;
    double labelwd = (o->figSizeX - hamwd)/2.0; /* Label width including gap. */
    double xpos = labelwd;
    double xneg = o->figSizeX - labelwd;

    /* Get number of levels {nneg,npos} in each tree, including root fraction: */
    double nneg = (dneg < 0 ? 0.0 : dneg + REL_NODE_WD);
    double npos = (dpos < 0 ? 0.0 : dpos + REL_NODE_WD);

    /* Divide plot width fairly among trees: */
    double xmid;
    if (nneg + npos == 0.0)
      { /* Nothing on both sides, {xmid} is irrelevant: */
        xmid = (xpos + xneg)/2;
      }
    else
      { double levelwd = hamwd/(nneg + npos);
        xmid = xpos + levelwd*npos;
      }
    (*xmidP) = xmid;
    (*xposP) = xpos;
    (*xnegP) = xneg;
  }

item_t *reverse_siblings(item_t *it)
  { item_t *r = NULL;
    while (it != NULL) { item_t *t = it->next; it->next = r; r = it; it = t; }
    return r;
  }

epswr_figure_t *init_plot_file
  ( char *dir, 
    char *prefix, 
    char *name, 
    int32_t seq, 
    char *suffix, 
    options_t *o
  )
  { double mm = epswr_pt_per_mm; /* One mm in pt. */
    double mrg = 4.0;  /* Default {epswr} margin width (pt). */
    double xtot = o->figSizeX*mm + 2.0*mrg;  /* Total figure width (pt). */
    double ytot = o->figSizeY*mm + 2.0*mrg;  /* Total figure height (pt). */
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( dir, prefix, name, seq, suffix, 
        xtot, ytot, mrg, mrg, mrg, mrg, eps_verbose
      );
    epswr_set_client_window(eps, 0, o->figSizeX, 0, o->figSizeY);
    epswr_set_label_font(eps, o->fontName, o->fontSize);
    return eps;
  }

void plot_slices
  ( epswr_figure_t *eps,  /* Postscript figure stream. */
    options_t *o,  /* Command line options. */
    double xbas,   /* X coordinate of uncut edge of ham. */
    double xtip,   /* X coordinate of cut edge of ham. */
    item_t *it,    /* The root node, or NULL. */
    int32_t dpmax,     /* Treat any nodes at this depth as leaves. */
    sign_t side    /* Side of plot. */
  )
  {
    /* Is there anything to plot? */
    if ((it == NULL) || (dpmax < 0)) { return; }
    
    /* Compute width {levelwd} of each level, and {nodewd} of each node: */
    double levelwd = fabs(xtip - xbas)/(dpmax + REL_NODE_WD);
    double nodewd = REL_NODE_WD * levelwd;

    /* Start of label: */
    double xlab = xtip - side*o->labelGap; /* Edge of label. */

    /* Pass 0 = paint interior, Pass 1 = draw outline, Pass 2 = write labels: */
    for (uint32_t pass = 0;  pass <= 2; pass++)
      {
        auto double do_plot ( item_t *itr, int32_t dp, double ybas );
          /* Plots the tree rooted at {itr}, assumed to be non-NULL and
            of depth {dp}. The bottom of the root and of the first leaf
            will be at ordinate {ybas}. Returns the top ordinate of the
            last leaf. */

        double do_plot ( item_t *itr, int32_t dp, double ybot )
          {
            assert(itr != NULL);
            assert(itr->val > 0);

            /* Plot node's rectangle: */
            bool_t is_leaf = ((dp == dpmax) || (itr->sub == NULL));
            double xini = xbas - side*dp*levelwd; /* Start X of node. */
            double xfin = (is_leaf ? xtip : xini - side*nodewd); /* Final X of node. */
            double ytop = ybot + itr->val/o->scale;

            /* Set corner Ys of slice: */
            if (pass == 0)
              { /* Fill slice: */
                epswr_rectangle(eps, xini,xfin, ybot,ytop, TRUE, FALSE);
              }
            else if (pass == 1)
              { /* Draw outline of slice: */
                epswr_segment(eps, xini, ybot, xfin, ybot);
                if (is_leaf || o->showNodes) { epswr_segment(eps, xfin, ybot, xfin, ytop); }
                if (o->showNodes && (dp > 0)) { epswr_segment(eps, xini, ybot, xini, ytop); }
                epswr_segment(eps, xfin, ytop, xini, ytop);
              }
            else if (pass == 2)
              { /* Do nothing */ }

            if (is_leaf)
              {
                /* Complete leaf node: */
                if (pass == 0)
                  { /* Do nothing: */ }
                else if (pass == 1)
                  { /* Do nothing: */ }
                else if (pass == 2)
                  { /* Write label and dot: */
                    double ylab = (ybot + ytop)/2;
                    epswr_dot(eps, xfin, ylab, 0.3, TRUE, FALSE);
                    double xalign = (side > 0 ? 1.0 : 0.0);
                    if ((o->fmt != NULL) && (strlen(o->fmt) > 0))
                      { char *tval = NULL;
                        char *tval = jsprintf(o->fmt, itr->val);
                        char *text = NULL;
                        if (side > 0)
                          { char *text = jsprintf("%s %s", itr->name, tval); }
                        else
                          { char *text = jsprintf("%s %s", tval, itr->name); }
                        epswr_label(eps, text, text, xlab, ylab, 0.0, TRUE, xalign, 0.5, TRUE, FALSE);
                        free(tval); free(text);
                      }
                    else
                      { epswr_label(eps, itr->name, itr->name, xlab, ylab, 0.0, TRUE, xalign, 0.5, TRUE,FALSE); }
                  }
                return ytop + o->leafGap;
              }
            else
              { /* Plot children branches: */

                double xjoin = xfin;                /* X of joined (parent) side of branch. */
                double xopen = xini - side*levelwd; /* X of open (child) side of branch. */

                double yjoin = ybot;  /* Bottom Y of branch, joined (parent) side. */
                double yopen = ybot;  /* Bottom Y of branch, open (child) side. */
                
                item_t *itp = itr->sub;
                while (itp != NULL)
                  {
                    double dy = itp->val/o->scale; /* Thickness of current slice. */

                    plot_branch_curved(eps, o, xjoin, xopen, yjoin, yopen, dy, pass);

                    /* Plot child and its subtree, get bottom Y of next slice at open end: */
                    yopen = do_plot(itp, dp+1, yopen);

                    /* Update bottom Y of next slice at joined end: */
                    yjoin += dy;

                    /* Get next item: */
                    itp = itp->next;
                  }
                return yopen;
              }
          }

        /* Set colors: */
        if (pass == 0)
          { float *c = o->fillColor.c;
            if (isnan(c[0]) || (c[0] == INF)) { /* Skip pass: */ continue; }
            epswr_set_fill_color(eps, c[0], c[1], c[2]);
          }
        else if (pass == 1)
          { float *c = o->drawColor.c;
            if (isnan(c[0]) || (c[0] == INF)) { /* Skip pass: */ continue; }
            epswr_set_pen(eps, c[0], c[1], c[2], 0.2, 0.0, 0.0);
          }
        else if (pass == 2)
          { float *c = o->textColor.c;
            if (isnan(c[0]) || (c[0] == INF)) { /* Skip pass: */ continue; }
            epswr_set_pen(eps, c[0], c[1], c[2], 0.2, 0.0, 0.0);
            epswr_set_fill_color(eps, c[0], c[1], c[2]);
          }

        /* Plot whole diagram: */
        double ymin = o->leafGap/2; /* Bottom Y of slices and ham. */
        do_plot(it, 0, ymin);
      }
  }

void plot_branch_polygonal
  ( epswr_figure_t *eps,  /* Postscript figure stream. */
    options_t *o,  /* Command line options. */
    double xjoin,  /* X coordinate of joined (parent) side of branch. */
    double xopen,  /* X coordinate of open (child) side of branch. */
    double yjoin,  /* Bottom Y of branch, joined (parent) side. */
    double yopen,  /* Bottom Y of branch, open (child) side. */
    double dy,     /* Y width of branch. */
    int32_t pass       /* Plotting pass. */
  )
  {
    double x[4], y[4]; /* Corners of connecting parallelogram. */

    /* Set abscissas of parallelogram: */
    x[0] = xjoin; x[1] = xopen; x[2] = xopen; x[3] = xjoin;

    /* Set corner Ys of slice: */
    y[0] = yjoin; y[1] = yopen; y[2] = yopen+dy; y[3] = yjoin+dy;

    if (pass == 0)
      { /* Fill slice: */
        epswr_polygon(eps, TRUE, x, y, 4, TRUE, FALSE, FALSE);
      }
    else  if (pass == 1)
      { /* Draw outline of slice: */
        epswr_segment(eps, x[0], y[0], x[1], y[1]);
        epswr_segment(eps, x[2], y[2], x[3], y[3]);
      }
    else if (pass == 2)
      { /* Do nothing: */ }
  }

void plot_branch_curved
  ( epswr_figure_t *eps,  /* Postscript figure stream. */
    options_t *o,  /* Command line options. */
    double xjoin,  /* X coordinate of joined (parent) side of branch. */
    double xopen,  /* X coordinate of open (child) side of branch. */
    double yjoin,  /* Bottom Y of branch, joined (parent) side. */
    double yopen,  /* Bottom Y of branch, open (child) side. */
    double dy,     /* Y width of branch. */
    int32_t pass       /* Plotting pass. */
  )
  {
    double x[8], y[8]; /* Bézier control points of lower and upper curves. */

    /* Set coordinates of corner points (ccw from lower join side): */
    x[0] = xjoin; y[0] = yjoin; 
    x[3] = xopen; y[3] = yopen; 
    x[4] = xopen; y[4] = yopen+dy;
    x[7] = xjoin; y[7] = yjoin+dy;

    /* Compute intermediate control points so that branch is horiz at corners: */
    for (int32_t i = 0; i < 8; i += 4)
      { x[i+1] = (2*x[i] + x[i+3])/3; y[i+1] = y[i];
        x[i+2] = (x[i] + 2*x[i+3])/3; y[i+2] = y[i+3];
      }

    if (pass == 0)
      { /* Fill slice: */
        epswr_bezier_polygon(eps, TRUE, x, y, 2, TRUE, FALSE, FALSE);
      }
    else  if (pass == 1)
      { /* Draw outline of slice: */
        epswr_curve(eps, x[0], y[0], x[1], y[1], x[2], y[2], x[3], y[3]);
        epswr_curve(eps, x[4], y[4], x[5], y[5], x[6], y[6], x[7], y[7]);
      }
    else if (pass == 2)
      { /* Do nothing: */ }
  }


void fini_plot_file(epswr_figure_t *eps, options_t *o)
  {
    epswr_end_figure(eps);
  }

options_t *get_options(int32_t argc, char **argv)
  {
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem");

    /* Defaults: */
    o->depth_pos = MAXDEPTH;
    o->depth_neg = MAXDEPTH;
    o->simplify = FALSE;
    o->figSizeX = DEFAULT_figSizeX;
    o->figSizeY = DEFAULT_figSizeY;
    o->plotWidth = -1.0; /* Default depends on {o->figSizeX}. */
    o->leafGap = DEFAULT_leafGap;
    o->showNodes = DEFAULT_showNodes;
    o->fillColor = DEFAULT_fillColor;
    o->drawColor = DEFAULT_drawColor;
    o->fontName = DEFAULT_fontName;
    o->fontSize = DEFAULT_fontSize;
    o->textColor = DEFAULT_textColor;
    o->fmt = DEFAULT_fmt;
    o->labelGap = DEFAULT_labelGap;
    o->scale = 0.0; /* Will be set automatically. */
    o->debug = FALSE;

    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Parse the options: */

    if (argparser_keyword_present(pp, "-depth"))
      { o->depth_pos = parse_depth(pp);
        o->depth_neg = parse_depth(pp);
      }

    o->simplify = argparser_keyword_present(pp, "-simplify");

    o->showNodes = argparser_keyword_present(pp, "-showNodes");

    o->debug = argparser_keyword_present(pp, "-debug");

    if (argparser_keyword_present(pp, "-figSize"))
      { o->figSizeX = argparser_get_next_double(pp, 1.0, 1000.0);
        o->figSizeY = argparser_get_next_double(pp, 1.0, 1000.0);
      }

    if (argparser_keyword_present(pp, "-plotWidth"))
      { o->plotWidth = argparser_get_next_double(pp, 0.0, 1000.0); }

    if (argparser_keyword_present(pp, "-leafGap"))
      { o->leafGap = argparser_get_next_double(pp, 0.0, 1000.0); }

    if (argparser_keyword_present(pp, "-fillColor"))
      { o->fillColor = parse_frgb(pp); }

    if (argparser_keyword_present(pp, "-drawColor"))
      { o->drawColor = parse_frgb(pp); }

    if (argparser_keyword_present(pp, "-font"))
      { o->fontName = argparser_get_next(pp);
        o->fontSize = argparser_get_next_double(pp, 0.0, 1000.0);
      }

    if (argparser_keyword_present(pp, "-textColor"))
      { o->textColor = parse_frgb(pp); }

    if (argparser_keyword_present(pp, "-fmt"))
      { o->fmt = argparser_get_next(pp); }

    if (argparser_keyword_present(pp, "-labelGap"))
      { o->labelGap = argparser_get_next_double(pp, 0.0, 1000.0); }

    if (argparser_keyword_present(pp, "-scale"))
      { o->scale = argparser_get_next_double(pp, 0.0, DBL_MAX); }

    argparser_finish(pp);

    /* Complete defaults: */
    if (o->plotWidth <= 0.0)
      { o->plotWidth = DEFAULT_plotWidthFraction * o->figSizeX; }

    return o;
  }

int32_t parse_depth(argparser_t *pp)
  {
    if (argparser_keyword_present_next(pp, "all"))
      { return MAXDEPTH; }
    else
      { return (int32_t)argparser_get_next_int(pp, -1, MAXDEPTH); }
  }

frgb_t parse_frgb ( argparser_t *pp )
  { frgb_t color;
    bool_t invisible = FALSE;
    for (uint32_t i = 0;  i < 3; i++)
      { color.c[i] = (float)argparser_get_next_double(pp, -1.0, +1.0);
        if (color.c[i] < 0.0) { invisible = TRUE; }
      }
    if (invisible)
      { return frgb_NoColor; }
    else
      { return color; }
  }

void data_error(int32_t line, char *msg)
  {
    fprintf(stderr, "%s:%d: **%s\n", "-", line, msg);
    exit(1);
  }
