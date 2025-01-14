#define PROG_NAME "make_test_classif_data"
#define PROG_DESC "make test data for vector classifiers"
#define PROG_VERS "1.0"

#define make_test_classif_data_C_COPYRIGHT \
  "Copyright � 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-21 11:55:46 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -problem {PROBLEM} \\\n" \
  "    -seed {SEED} \\\n" \
  "    -samples {NS} \\\n" \
  "    -prefix {PREFIX} \\\n" \
  "    [ -attributes {NAR} {NAI} ] \\\n" \
  "    [ -classes {NC} ] \\\n" \
  "    [ -noise {SIGMA} ] \\\n" \
  "    [ -verify { N | Y } ] \\\n" \
  "    [ -image {IMG_SIZE} {SUBSMP} ] \\\n" \
  "    [ -grid { N | Y } ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Writes test data files for vector classifiers, containing {NS} samples" \
  " from a classification problem specified by {PROBLEM}, {NAR} and {NC}.\n" \
  "\n" \
  "  Each sample is a vector of attributes.  Each attribute is a real number.  The total number of" \
  " attributes per sample is {NA = NAR + NAI} where {NAR} and {NAI} are" \
  " user-specifiable numbers.  The first {NAR} attributes are \"relevant\" to the test problem;" \
  " the the remaining {NAI} attributes are \"irrelevant\" random numbers, uniformly distributed" \
  " in {U = [-1 _ +1]}, appended to them.\n" \
  "\n" \
  "   Considering only the {NAR} relevant attributes, the space of all valid" \
  " attribute vectors consists of {NC}" \
  " disjoint /class domains/ contained in {U^NAR}, numbered {1..NC}.  The program" \
  " generates {NS} samples in these domains, and assigns to each sample" \
  " a /nominal class/ which is the index of the containing domain.  The samples" \
  " may be generated at random or may be taken from a regular grid. The number" \
  " and shape of the domains, and the random sampling probability distribution" \
  " in each domain depends on the {PROBLEM} and" \
  " on the user-specified {NAR} and {NC} parameters.  For some {PROBLEM}s" \
  " the number of relevant attributes {NAR} and/or classes {NC} may be restricted or fixed.  The class domains" \
  " may not cover the whole of {U^NAR} and may or may not touch each other. Some" \
  " domains may be zero-measure sets with dimension smaller than {NAR}.\n" \
  "\n" \
  "  If grid sampling is requested, the number {NS} is implicitly rounded up to the next" \
  " perfect power {NG^NAR}, and the problem is sampled at a regular grid of points" \
  " with {NAR} samples along each axis and spanning the hypercube {U^NAR}.  Only" \
  " those samples that fall inside one of" \
  " the class domains will be written out.  In that case the actual number of" \
  " samples in the output may be substantially larger or smaller than {NS}.\n" \
  "\n" \
  "  After generating each sample (whether randomly or from a grid) and" \
  " appending the {NAI} irrelevant attributes, the program" \
  " adds to each attribute (relevant or irrelevant)" \
  " a Gaussian noise with mean 0 and deviation {SIGMA}.  The noise is truncated to {4*SIGMA}" \
  " so that each final attribute is in the range {V = [-VMAX _ +VMAX]}" \
  " where {VMAX = 1 + 4*SIGMA}.  The nominal class of the sample is retained," \
  " even though the perturbation may cause the" \
  " attribute vector to fall outside the domain of that class and" \
  " possibly into the domain of a different class.\n" \
  "\n" \
  "OUTPUTS\n" \
  "  The program writes the samples generated in each class" \
  " to a file \"{PREFIX}-c{CCC}.dat\" where {CCC} is the class" \
  " index.  The file format is defined by {rn_classif_dataset_write}.\n" \
  "\n" \
  "  The program may also write a PPM image file showing the\n" \
  " classes.  See the \"-image\" option for details.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -problem {PROBLEM}\n" \
  "    This mandatory argument specifies the ideal class domains and" \
  " the sampling probability distributions.\n" \
  "\n" \
  "  -samples {NS}\n" \
  "    This mandatory argument specifies the number of samples to generate.\n" \
  "\n" \
  "  -seed {SEED}\n" \
  "    This mandatory argument specifies a seed for the random number generator.\n" \
  "\n" \
  "  -prefix {PREFIX}\n" \
  "    This mandatory argument specifies the prefix for all output filenames.\n" \
  "\n" \
  "  -attributes {NAR} {NAI}\n" \
  "    This optional argument specifies the number of attributes in" \
  " each sample, namely {NAR} relevant ones and {NAI} irrelevant ones.  The allowed" \
  "  values of {NAR} depend on the {PROBLEM}.  If {NAR} is zero" \
  " (the default) the number of relevant attributes is selected by the {PROBLEM}.  If {NAI} is zero" \
  " (the default) there are no irrelevant attributes.\n" \
  "\n" \
  "  -classes {NC}\n" \
  "    This optional argument specifies the number of classes into which" \
  " the space of all attribute vectors is partitioned.  If {NC} is zero" \
  " (the default) the number of classes is selected by the {PROBLEM}.  The" \
  " allowed range of {NC} depends on the {PROBLEM}.\n" \
  "\n" \
  "  -noise {SIGMA}\n" \
  "    This optional argument specifies the standard deviation" \
  " of the Gaussian perturbation to be added to each attribute" \
  " after sampling.  Note that the noise may move some samples" \
  " into the domain of a different class.  The default is {SIGMA=0}" \
  " (no perturbation).\n" \
  "\n" \
  "  -verify { N | Y }\n" \
  "    This optional argument specifies whether the procedure should" \
  " check that each randomly generated sample belongs to the domain of its nominal" \
  " class.  It has no effect on grid sampling.  Errors before the" \
  " random perturbation are fatal; errors after the perturbation are" \
  " counted and reported to {stderr}.  The default is \"N\" (no checking).\n" \
  "\n" \
  "  -grid { N | Y }\n" \
  "    This optional argument specifies whether the procedure should" \
  " generate sample the domains at a regular grid instead of" \
  " randomly.  The default is \"N\" (random sampling).\n" \
  "\n" \
  "  -image {SIZE} {SUBSMP}\n" \
  "    This optional argument directs the program to write a PPM" \
  " image file called \"{PREFIX}.ppm\" showing the domain" \
  " classes.  This option can be used" \
  " only when {NAR=2}. The image will have {SIZE} columns and {SIZE}" \
  " rows and will span the augmented square {V^2}.  Each domain will" \
  " be painted with a distinct nonwhite color, with brightness" \
  " increasing from 1 to {NC}.  The complement of all domains" \
  " is painted white.  The program will generate {SUBSMP^2} samples" \
  " inside each pixel, and average their colors to obtain the pixel" \
  " color.  If \"-noise\" is requested, each sample is perturbed by" \
  " the noise, and then the original sample is painted with the color" \
  " of the domain that contains the perturbed point.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  The San Francisco Exploratorium.\n" \
  "\n" \
  "SMELL ALSO\n" \
  "  The roses along the road.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2010-05-21 by J. Stolfi: created program.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " make_test_classif_data_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <uint16_image.h>
#include <float_image.h>
#include <r2.h>
#include <rn.h>
#include <rn_classif.h>
#include <rn_classif_image.h>
#include <rn_classif_test.h>
  
#define MAX_ATTRIBS 128
  /* Max num of relevant attributes (param safety only) */

#define MAX_CLASSES 128
  /* Max num of classes (param safety only) */

#define MAX_SAMPLES 100000000
  /* Max num of samples to generate (param safety only) */

#define MAX_SEED (~0LU)
  /* Max seed value (param safety only) */

#define MAX_NOISE 10.0
  /* Max variance of post-sampling noise (param safety only) */

#define MAX_IMAGE_SIZE 1024
  /* Max rows and columns in image (param safety only) */

#define MAX_SUBSMP 5
  /* Max pixel subsampling order for image (param safety only) */

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { char *problem;      /* Name of classification problem. */
    int attributes_rel; /* Number of relevant attributes desired. */
    int attributes_irr; /* Number of irrelevant attributes to append. */
    int classes;        /* Number of classes desired. */
    int samples;        /* Number of samples to generate. */
    uint32_t seed;      /* Seed for randomizer. */
    double noise;       /* Deviation of post-sampling noise. */
    bool_t verify;      /* TRUE to verify the generated samples. */
    bool_t grid;        /* TRUE to sample points on a regular grid. */
    char *prefix;      /* The ouput file prefix. */
    /* Image output options: */
    int image_size;     /* Height and width of PPM image, or 0 if none. */
    int image_subsmp;   /* Will generate {subsmp*subsmp} samples in each pixel. */
  } options_t;

/* PROBLEMS */

typedef enum
  { problem_kind_SATURN, /* The "saturn" ("B9") dataset from Papa et al.(2009); {NC = 2}. */
    problem_kind_PETALS, /* The "petals" ("B10") dataset from Papa et al.(2009); {NC = 4}. */
    problem_kind_VESSEL, /* The "boat" ("B11") dataset from Papa et al.(2009); {NC = 3}. */
    problem_kind_MBALLS, /* Each class {2..NC} comprises {(NC-1)^(NAR-1)} balls; class 1 is background. */
    problem_kind_SHELLS, /* Classes {2..NC} are concentric shells; class 1 is background and center. */
    problem_kind_NUMBER  /* Number of valid problem kinds. Must be last. */
 } problem_kind_t;
 /* Numeric problem kinds.  The valid kinds are {0..problem_kind_NUMBER-1} */

char *problem_kind_name[problem_kind_NUMBER] = 
  { [problem_kind_SATURN] = "saturn", 
    [problem_kind_PETALS] = "petals", 
    [problem_kind_VESSEL] = "vessel", 
    [problem_kind_MBALLS] = "mballs", 
    [problem_kind_SHELLS] = "shells", 
  };
  /* External names of the problem kinds. */
  
typedef void problem_def_proc_t(int *NAR, int *NC);
  /* Type of a problem definition procedure.  It sets {*NC}
    and {*NAR} to suitable defaults if any of them are zero, then 
    performs range checking. */

/* INTERNAL PROTOTYPES */

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);

rn_classif_problem_t *get_problem(char *name, int NAR, int NC);
  /* Returns a {rn_classif_problem_t} given the problems {name} and the desired
    numbers of attributes {NAR} and classes {NC}.  Fails if {NAR,NC}
    are not a valid combination for that problem. */
    
void generate_dataset(rn_classif_problem_t *P, options_t *o, rn_classif_dataset_t **DP, int **classDP);
  /* Returns a {dataset_t} containing samples of the classification
    problem {P} according to options {o}. */

void generate_raw_dataset_random(rn_classif_problem_t *P, int NA, int NS, uint32_t seed, bool_t verify, rn_classif_dataset_t **DP, int **classDP);
void generate_raw_dataset_grid(rn_classif_problem_t *P, int NA, int NS, uint32_t seed, rn_classif_dataset_t **DP, int **classDP);
  /* Returns a {dataset_t} containing random or grid samples of the classification
    problem {P} according to options {o}.  Only the {P->NAR} relevant attributes
    of each sample are set.  Does not add any noise. */
   
void dataset_stats_print(FILE *wr, int NS, int NC, int classD[]);
  /* Prints statistics on number and percentage of samples per 
    class in the dataset {D}. */

void output_dataset(rn_classif_dataset_t *D, int classD[], int NC, options_t *o);
  /* Writes the dataset as {NC} files "{o->prefix}-c{NNN}.dat" where
    {CCC} is each class in 3 digit format. See {rn_classif_dataset_write}
    for the format. */
    
void verify_dataset(rn_classif_problem_t *P, rn_classif_dataset_t *D, int classD[]);
  /* Compares the classification {classD} for the samples of {D}
    with the classification specified by {P.lab} with {P.NA} atributes
    and {P.NC} classes.  Assumes that {classD} has only classes in {0..P.NC} too. */

void output_image(rn_classif_problem_t *P, options_t *o);
  /* Writes a PPM file "{o->prefix}.ppm" showing the class domains of problem {P}. */

void add_random_noise(rn_classif_dataset_t *D, double sigma, uint32_t seed);
  /* Add the post-sampling noise.  The perturbations are a fixed 
    random function of the seed, indepednent of the random
    values used in sampling. */

void append_irrelevant_attributes(rn_classif_dataset_t *D, rn_classif_problem_t *P, uint32_t seed);
  /* Appends {D.NA - P.NAR} random values uniform in {U} to each 
    sample in {D}.  The values are a fixed 
    random function of the seed, independent of those used in 
    sampling and noise. */

void get_grid_indices(int g, int NG, int NAR, int ix[]);
  /* The integer {g} must be in the range {0..NG^NAR-1}. Decomposes
     {g} into its {NAR} digits {ix[0..NAR-1]} in base {NG}. */

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
    
    rn_classif_problem_t *P = get_problem(o->problem, o->attributes_rel, o->classes);
    
    rn_classif_dataset_t *D;
    int *classD;
    generate_dataset(P, o, &D, &classD);
    
    if (o->verify) { verify_dataset(P, D, classD); }
    
    output_dataset(D, classD, P->NC, o);
    
    if (o->image_size > 0) { output_image(P, o); }
    
    return 0;
  }

void generate_dataset(rn_classif_problem_t *P, options_t *o, rn_classif_dataset_t **DP, int **classDP)
  { 
    int NA = P->NA + o->attributes_irr;
    if (o->grid)
      { generate_raw_dataset_grid(P, NA, o->samples, o->seed, DP, classDP); }
    else
      { generate_raw_dataset_random(P, NA, o->samples, o->seed, o->verify, DP, classDP); }
    rn_classif_dataset_t *D = (*DP);   
    int *classD = (*classDP);
    dataset_stats_print(stderr, D->NS, P->NC, classD);
    append_irrelevant_attributes(D, P, o->seed);
    if (o->noise > 0) { add_random_noise(D, o->noise, o->seed); } 
  }

void dataset_stats_print(FILE *wr, int NS, int NC, int classD[])
  {
    int *num = notnull(malloc((NC+1)*sizeof(int)), "no mem");
    rn_classif_class_count(NS, classD, NC, num); 
    rn_classif_class_count_print(wr, NC, num);
    free(num);
  }

void output_dataset(rn_classif_dataset_t *D, int classD[], int NC, options_t *o)
  { int cl;
    for (cl = 1; cl <= NC; cl++)
      { char *fname = NULL;
        char *fname = jsprintf("%s-c%03d.dat", o->prefix, cl);
        FILE *wr = open_write(fname, TRUE);
        rn_classif_dataset_write(wr, D, cl, classD);
        fclose(wr);
        free(fname);
      }
  }

void verify_dataset(rn_classif_problem_t *P, rn_classif_dataset_t *D, int classD[])
  { int *classX = notnull(malloc(D->NS*sizeof(int)), "no mem");
    rn_classif_dataset_label(D, P->NA, P->NC, P->lab, classX);
    int *ncc;
    rn_classif_cross_matrix_build(D->NS, classD, P->NC, classX, P->NC, &ncc);
    rn_classif_cross_matrix_print(stderr, P->NC, P->NC, ncc, TRUE);
    free(ncc);
  }

void output_image(rn_classif_problem_t *P, options_t *o)
  { r2_t ctr = (r2_t){{ 0,0, }};   /* Center of imaged area. */
    double HV = 1.0 + 4*o->noise;  /* Half-extent of imaged area. */
    srandom(o->seed + 314159);
    uint16_image_t *img = rn_classif_compute_image(P->NA, P->NC, P->lab, &ctr, HV, o->image_size, o->image_subsmp, o->noise);
    char *fname = jsprintf("%s.ppm", o->prefix);
    uint16_image_write_pnm_named(fname, img, FALSE, TRUE);
    uint16_image_free(img);
    free(fname);
  }
    
void append_irrelevant_attributes(rn_classif_dataset_t *D, rn_classif_problem_t *P, uint32_t seed)
  {
    int NS = D->NS;
    int NAR = P->NA;
    int NA = D->NA;
    /* Append the irrelevant attributes {NAR..NA-1}: */
    /* Must do this separately so that the essential ones are fixed for a given {seed}. */
    srandom(seed + 4615);
    int i, t;
    for (i = 0; i < NS; i++)
      { double *attri = D->smp[i];
        for (t = NAR; t < NA; t++) { attri[t] = 2*drandom() - 1; }
      }
  }

void add_random_noise(rn_classif_dataset_t *D, double sigma, uint32_t seed)
  {
    int NS = D->NS;
    int NA = D->NA;
    srandom(seed + 19501129);
    int i, t;
    for (i = 0; i < NS; i++)
      { double *attri = D->smp[i];
        for (t = 0; t < NA; t++) 
          { /* Add to {attri[t]} a Gaussian deviate truncated to {4*sigma} */
            double dit = sigma*fmax(-4.0, fmin(+4.0, dgaussrand()));
            attri[t] += dit;
          }
      }
  }
     
void generate_raw_dataset_random(rn_classif_problem_t *P, int NA, int NS, uint32_t seed, bool_t verify, rn_classif_dataset_t **DP, int **classDP)
  {
    int NC = P->NC;
    int NAR = P->NA;
    rn_classif_dataset_t *D = rn_classif_dataset_new(NS, NA);
    int *classD = notnull(malloc(NS*sizeof(int)), "no mem");

    /*  Generate samples: */
    srandom(seed);
    int nerr = 0;
    int i;
    for (i = 0; i < NS; i++)
      { double *attri = notnull(malloc(NA*sizeof(double)), "no mem");
        D->smp[i] = attri;
        P->gen(i, NAR, NC, attri, &(classD[i]));
        if (verify) 
          { int ver = P->lab(NAR, NC, attri);
            if (classD[i] != ver)
              { /* Since sample is unperturbed,this is a fatal error or very bad luck: */
                fprintf(stderr, "** inconsistent class %d != %d for sample %d:\n", classD[i], ver, i);
                rn_classif_sample_print(stderr, "  ", NAR, attri, " ", "\n");
                nerr++;
                if (nerr > 100) { fprintf(stderr, "** too many errors\n"); exit(1); }
              }
          }
      }
    (*DP) = D;
    (*classDP) = classD;
  }
  
void generate_raw_dataset_grid(rn_classif_problem_t *P, int NA, int NS, uint32_t seed, rn_classif_dataset_t **DP, int **classP)
  {
    int NC = P->NC;
    int NAR = P->NA;
    
    /* Compute the grid order {NG} such that {NG^NAR >= NS}: */
    int NG = (int)ceil(pow(NS, 1.0/NAR));
    int NT = ipow(NG, NAR); /* Total number of tentative samples. */
    assert(NT >= NS);
    demand(NT <= MAX_SAMPLES, "too many samples in grid");
    
    /* Allocate sample array with maximum size: */
    int *classD = notnull(malloc(NT*sizeof(int)), "no mem");
    rn_classif_dataset_t *D = rn_classif_dataset_new(NT, NA);
    
    /* Each sample is identified by a tuple of {NAR} indices in {0..NG-1}: */
    int ix[NAR];
    
    /* Generate grid samples and store those that have nonzero class. */
    /* Those samples are stored in postions {0..NOK-1} of {attr} and {classD}. */
    NS = 0;
    srandom(seed);
    double attr[NA]; /* Temporary attribute vector. */
    int g;
    for (g = 0; g < NT; g++)
      { /* Break {g} down into {NAR} grid indices: */
        get_grid_indices(g, NG, NAR, ix);
        /* Convert the indices into attributes in {U}: */
        int t;
        for (t = 0; t < NAR; t++) { attr[t] = 2*(ix[t] + 0.5)/NG - 1; }
        /* Now use the problem's classifier to obtain the classD: */
        int cl = P->lab(NAR, NC, attr);
        if (cl != 0)
          { /* Sample is inside some domain, store it: */
            double *attrn = notnull(malloc(NA*sizeof(double)), "no mem");
            for (t = 0; t < NAR; t++) { attrn[t] = attr[t]; }
            D->smp[NS] = attrn;
            classD[NS] = cl;
            NS++;
          }
      }
      
    /* Trim excess storage and return: */
    D->NS = NS;
    D->smp = notnull(realloc(D->smp, NS*sizeof(double*)), "no mem");
    classD = notnull(realloc(classD, NS*sizeof(int)), "no mem");
    (*DP) = D;
    (*classP) = classD;
  }
  
void get_grid_indices(int g, int NG, int NAR, int ix[])
  {
    int k;
    int tmp = g;
    for (k = 0; k < NAR; k++) { ix[k] = tmp % NG; tmp = tmp/NG; } 
  }
  
rn_classif_problem_t *get_problem(char *name, int NAR, int NC)
  {
    problem_kind_t kind;
    for (kind = 0; kind < problem_kind_NUMBER; kind++) 
      { if (strcmp(name, problem_kind_name[kind]) == 0)
          { /* Found the problem: */
            /* Get the problem-specific procedure {def}: */
            problem_def_proc_t *chk = NULL;
            rn_classif_labeler_t *lab = NULL;
            rn_classif_thrower_t *gen = NULL;
            switch(kind)
              {
                case problem_kind_SATURN: 
                  chk = &rn_classif_test_check_saturn;
                  lab = &rn_classif_test_label_saturn;
                  gen = &rn_classif_test_throw_saturn;
                  break;
                case problem_kind_PETALS: 
                  chk = &rn_classif_test_check_petals;
                  lab = &rn_classif_test_label_petals;
                  gen = &rn_classif_test_throw_petals;
                  break;
                case problem_kind_VESSEL: 
                  chk = &rn_classif_test_check_vessel;
                  lab = &rn_classif_test_label_vessel;
                  gen = &rn_classif_test_throw_vessel;
                  break;
                case problem_kind_MBALLS: 
                  chk = &rn_classif_test_check_mballs;
                  lab = &rn_classif_test_label_mballs;
                  gen = &rn_classif_test_throw_mballs;
                  break;
                case problem_kind_SHELLS: 
                  chk = &rn_classif_test_check_shells;
                  lab = &rn_classif_test_label_shells;
                  gen = &rn_classif_test_throw_shells;
                  break;
                default: assert(FALSE);
              }
            rn_classif_problem_t *P = notnull(malloc(sizeof(rn_classif_problem_t)), "no mem");
            (*P) = (rn_classif_problem_t) { .NA = NAR, .NC = NC, .lab = lab, .gen = gen };
            chk(&(P->NA), &(P->NC));
            return P;
          }
      }
    demand(FALSE, "unrecognized problem kind");
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
    
    argparser_get_keyword(pp, "-problem");
    o->problem = argparser_get_next(pp);

    argparser_get_keyword(pp, "-samples");
    o->samples = argparser_get_next_int(pp, 1, MAX_SAMPLES);

    argparser_get_keyword(pp, "-seed");
    o->seed = argparser_get_next_uint(pp, 1, MAX_SEED);
    
    argparser_get_keyword(pp, "-prefix");
    o->prefix = argparser_get_next(pp);
    
    if (argparser_keyword_present(pp, "-attributes"))
      { o->attributes_rel = argparser_get_next_int(pp, 0, MAX_ATTRIBS);
        o->attributes_irr = argparser_get_next_int(pp, 0, MAX_ATTRIBS - o->attributes_rel); 
      }
    else
      { o->attributes_rel = 0; o->attributes_irr = 0; }

    if (argparser_keyword_present(pp, "-classes"))
      { o->classes = argparser_get_next_int(pp, 0, MAX_CLASSES); }
    else
      { o->classes = 0; }

    if (argparser_keyword_present(pp, "-noise"))
      { o->noise = argparser_get_next_double(pp, 0, MAX_NOISE); }
    else
      { o->noise = 0; }

    if (argparser_keyword_present(pp, "-verify"))
      { o->verify = argparser_get_next_bool(pp); }
    else
      { o->verify = FALSE; }

    if (argparser_keyword_present(pp, "-grid"))
      { o->grid = argparser_get_next_bool(pp); }
    else
      { o->grid = FALSE; }

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
