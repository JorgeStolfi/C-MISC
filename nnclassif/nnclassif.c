#define PROG_NAME "nnclassif"
#define PROG_DESC "nearest-neighbor classifier with opf and other metrics"
#define PROG_VERS "1.0"

#define nnclssif_C_COPYRIGHT \
  "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-21 11:55:26 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -model {MODEL_FILE}" \
  "    -check {CHECK_FILE}" \
  "    [ -train {TRAIN_FILE}" ] \
  "    -seed {SEED} \\\n" \
  "    -prefix {PREFIX} \\\n" \
  "    [ -opf { N | Y } ] \\\n" \
  "    [ -image {IMG_SIZE} ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program bla bla bla" \
  " bla bla bla bla {X+Y} bla bla" \
  " bla {INFILE} bla \"foobar.ppm\" bla bla bla\n" \
  "\n" \
  "  Beware that bla bla bla BLEBBLE BLOB bla" \
  " bla bla bla bla.\n" \
  "\n" \
  "OPTIONS\n" \
  "  \"Freedom does not mean having options, it means" \
  " not having to opt.\"\n" \
  "\n" \
  "  -model {MODEL_FILE}\n" \
  "    This mandatory argument specifies the data file that contains" \
  " the representative samples used by the classifier.\n" \
  "\n" \
  "  -check {CHECK_FILE}\n" \
  "    This mandatory argument specifies the data file that contains" \
  " the samples to be classified\n" \
  "\n" \
  "  -train {TRAIN_FILE}\n" \
  "    This optional argument specifies the data file that contains" \
  " additional samples that can be used to improve the classifier" \
  " before checking.\n" \
  "\n" \
  "  -seed {SEED}\n" \
  "    This mandatory argument specifies a seed for the random number generator.\n" \
  "\n" \
  "  -prefix {PREFIX}\n" \
  "    This mandatory argument specifies the prefix for all output filenames.\n" \
  "\n" \
  "  -opf { N | Y }\n" \
  "    This optional argument requests the use of handicaps" \
  " when computing distances to the representantive samples" \
  " in the classifer.  The handicaps are" \
  " computed by the OPF algorithm   The default is \"N\" (no OPF handicaps).\n" \
  "\n" \
  "  -image {SIZE} {SUBSMP}\n" \
  "    This optional argument directs the program to write a PPM" \
  " image showing the domain classes defined by the classifier.  This" \
  " option can be used only when the number of attributes {NA} is 2. The" \
  " image will have {SIZE} columns and {SIZE}" \
  " rows and will span the square {V^2} where {V = [-1.2 _ +1.2]}.  Each domain will" \
  " be painted with a distinct nonwhite color, with brightness" \
  " increasing from 1 to {NC}.  The complement of all domains" \
  " is painted white.  The program will generate {SUBSMP^2} samples" \
  " inside each pixel, and average their colors to obtain the pixel" \
  " color.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  Farther than others, by standing on their shoulders.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2010-05-23 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " nnclassif_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <uint16_image.h>
#include <float_image.h>
#include <rn.h>
#include <rn_classif.h>
#include <jsrandom.h>
/* #include <rn_classif_image.h> */

#define MAX_SEED (~0LU)
  /* Max seed value (param safety only) */

#define MAX_IMAGE_SIZE 1024
  /* Max rows and columns in image (param safety only) */

#define MAX_SUBSMP 5
  /* Max pixel subsampling order for image (param safety only) */

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { char *model;          /* Name of test dataset. */
    char *train;          /* Name of test dataset. */
    char *check;          /* Name of test dataset. */
    uint32_t seed;        /* Seed for randomizer. */
    bool_t opf;           /* TRUE uses OPF handicaps. */
    char *prefix;      /* The ouput file prefix. */
    /* Image output options: */
    int image_size;     /* Height and width of PPM image, or 0 if none. */
    int image_subsmp;   /* Will generate {subsmp*subsmp} samples in each pixel. */
  } options_t;

/* INTERNAL PROTOTYPES */

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    rn_classif_dataset_t *M = rn_classif_dataset_read(o->model);
    int NA = M->NA; /* Number of attributes (including irrelevant ones). */
    int NC = M->NC; /* Number of classes. */
   
    rn_classif_dataset_t *C = rn_classif_dataset_read(o->check);
    demand(C->NA == NA, "wrong num attributes");
    demand(C->NC == NC, "wrong num classes");
   
    /* Allocate and compute the handicaps: */
    double *H = compute_handicaps(&M, o);
    
    /* LOCAL PROC PROTOTYPES */
    
    auto double hdist(int NA, double p[], rn_classif_dataset_t *M, int im);
      /* Computes the distance from {p[0..NA-1]} to
        representative {ip} of the model {M}. The 
        distace is modified by the handicap of {im}, if any. */
    
    auto void nnc(int i, int NA, double p[], int NC, int *class);
      /* The nearest neighbor classifiers with model {M} and
        handicaps {H}, packaged as a {rn_classif_proc_t}. */
     
    /* Train the model if so requested: */    
    if (o->train != NULL)
      { rn_classif_dataset_t *T = rn_classif_dataset_read(o->train);
        demand(T->NA == NA, "wrong num attributes");
        demand(T->NC == NC, "wrong num classes");
        
        train(&M, &T);
      }
    
    /* Evaluate the model on the dataset {C}: */    
    rn_classif_dataset_check(NA, nnc, &C);
    
    /* Write image of the {nn} classifier if so requested: */    
    if (o->image_size > 0) 
      { r2_t ctr = (r2_t){{ 0,0, }};   /* Center of imaged area. */
        double HV = 1.2;  /* Half-extent of imaged area. */
        uint16_image_t *img = rn_classif_compute_image(NA, nnc, NC, &ctr, HV, o->image_size, o->image_subsmp, 0.0);
        uint16_image_write_pnm_named(o->image_name, img, FALSE, TRUE);
        uint16_image_free(img);
      }
    
    return 0;
    
    /* LOCAL PROC IMPLEMENTATIONS */
    
    double hdist(int NA, double p[], int im)
      {
        assert((im >= 0) && (im < NS));
        double *q = &(M->attr[im*M->NA]);
        double edist = rn_dist(NA, p, q);
        return fmax(edist, H[im]);
      }

    void nnc(int i, int NA, double p[], int NC, int *class)
      { 
        assert(i < 0);
        rn_classif_nearest(&M, hdist, NA, p, NC, class)
      }
    
  }
 
!!!
  
void hproc(){
if (method == "1nn") {
      /* That is it. */
    } else if (method == "opf") {

   } else {
      affirm(0, "invalid method");
    }
}

/* Classes of {T} computed by {(M,H,classM)}: */
    int classX = notnull(malloc(T->NS*sizeof(int)), "no mem");
    
    maxIter = 3*nZ[2]
    free(fail);

void rn_classif_dataset_stats_print(FILE *wr, rn_classif_dataset_t *D);
  /* Prints statistics on number and percentage of samples per 
    class in the dataset {D}. */

void rn_classif_dataset_stats_print(FILE *wr, rn_classif_dataset_t *D)
  {
    int *num = notnull(malloc((D->NC+1)*sizeof(int)), "no mem");
    rn_classif_dataset_class_count(D, num); 
    rn_classif_dataset_class_count_print(wr, D->NC, num);
    free(num);
  }

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-model");
    o->model = argparser_get_next(pp);

    argparser_get_keyword(pp, "-check");
    o->check = argparser_get_next(pp);
    
    if (argparser_keyword_present(pp, "-train"))
      { o->train = argparser_get_next(pp); }
    else
      { o->train = NULL; }

    argparser_get_keyword(pp, "-seed");
    o->seed = argparser_get_next_uint(pp, 1, MAX_SEED);
    
    argparser_get_keyword(pp, "-prefix");
    o->prefix = argparser_get_next(pp);
    
    if (argparser_keyword_present(pp, "-opf"))
      { o->opf = argparser_get_next_bool(pp); }
    else
      { o->opf = FALSE; }

    if (argparser_keyword_present(pp, "-image"))
      { o->image_size = argparser_get_next_int(pp, 0, MAX_IMAGE_SIZE);
        o->image_subsmp = argparser_get_next_int(pp, 0, MAX_SUBSMP);
      }
    else
      { o->image_size = 0;
        o->image_subsmp = 0;
      }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }


//  /* GAWK
//  
//  
//  BEGIN {
//  
//    USAGE = "test-1nn.gawk -v dataset={NAME} -v method={METHOD}"
//    
//    /* where {METHOD} is "opf" or "1nn". */
//  
//    /* Global variables: */
//    /* {m} Number of classes; classes are numbered {1..m}. */
//    /* {NS} Total number of points. Points are numbered {0..NS-1}. */
//    /* {nC[1..m]} is the number of samples in each class. */
//    /* {X[0..NS-1]} X-coordinate of each sample. */
//    /* {Y[0..NS-1]} Y-coordinate of each sample. */
//    /* {class[0..NS-1]} class of each sample. */
//    /* {part[0..NS-1]} part of each sample: 1 for {Z1}, 2 for {Z2}, 3 for {Z3}. */
//    /* {fail[0..NS-1]} set to 1 if classification failed, 0 if suceeded, -1 if not tried. */
//    /* {H[0..NS-1]} Handicap for distance classification. */
//    /* {nZ[1..3]} is the number of samples in each part of the dataset. */
//    /* {repr[0..nZ[1]-1]} are the elements of part {Z1}. */
//    
//    for (adjust = 0; adjust <= 1; adjust++) {
//      ntries = 10;
//      rS_sum = 0;
//      for (try = 0; try < ntries; try++) {
//        verbose = (try == 0);
//        dump = (try == 0);
//        run_test(method,verbose,dump,adjust);
//        if (verbose) { fprintf(stderr, "classifying Z3 with Z1 ... " ); }
//        classify_points(3,method,verbose);
//        rS = nS/(nS + nF);
//        fprintf("trial %3d  success rate = %6.4f\n", try, 100*rS > "/dev/stderr"
//        rS_sum += rS);
//      }
//      rS = rS_sum/ntries;
//      xadj = (adjust ? "with" : "without" );
//      fprintf("mean success rate %s learning = %6.4f\n", xadj, 100*rS > "/dev/stderr"
//    }
//  }
//  
//  function run_test(method,verbose,dump,adjust,  cl,i,pt) {
//  
//    /* If {verbose} is true writes diagnostics */
//    /* if {dump} is true writes the datasets for plotting. */
//    /* If {adjust} tries to adjust {Z1} for best classif in {Z2}. */
//  
//    m = num_classes(dataset));  /* Number of classes. */
//    split("", nC);             /* {nC[cl]} is the number of samples in class {cl}. */
//  
//    for (cl = 1; cl <= m; cl++) { nC[cl] = 0; }
//  
//    /* Randomly partition points into {Z1,Z2,Z3}: */
//    for (i = 0; i < NS; i++) {
//      coin = myrand();
//      if (coin < 0.30) { 
//        pt = 1;
//      } else if (coin < 0.50) {
//        pt = 2;
//      } else {
//        pt = 3;
//      }
//      part[i] = pt;
//    }
//  
//   /* Classify the points of {Z3} using {Z1}: */
//    split("", fail);
//    if (verbose) { fprintf(stderr, "classifying Z3 with Z1 ... " ); }
//    classify_points(3,method,verbose);
//    if (dump) { write_datasets(method,"ini"); }
//     
//    if (adjust) {
//      /* Reclassify {Z3} using {z1}: */
//      adjust_representatives(method,verbose); 
//  
//      /* Reclassify {Z3} using {z1}: */
//      if (verbose) { fprintf(stderr, "\n" ); }
//      if (verbose) { fprintf(stderr, "classifying Z3 with Z1 ... " ); }
//      classify_points(3,method,verbose);
//      if (dump) { write_datasets(method,"opt"); }
//    }
//  }
//  
//  function dist(i1,i2,  dx,dy) {
//    /* Euclidean distance between points {i1} and {i2}. */
//    dx = X[i1] - X[i2];
//    dy = Y[i1] - Y[i2];
//    return sqrt(dx*dx + dy*dy);
//  }
//    
//  */
