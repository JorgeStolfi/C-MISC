
#define gdr_quantiles_HEADERS  "lo 95%", "median", "hi 95%" 

    char *tit_n[gdr_quantiles_NUM] = { gdr_quantiles_HEADERS };

    fprintf(stderr, "%6s  %32s  %32s", "", "population", "surviving lineages");
    fprintf(stderr, "\n");

    char *dash32 = "--------------------------------";
    fprintf(stderr, "%6s  %32s  %32s", "", dah32, dash32);
    fprintf(stderr, "\n");

    fprintf(stderr, "%6s", "year");
    for (int32_t m = 0; m < 2; m++)
      { fprintf(stderr, " ");
        for (int32_t q = 0; q < nq; q++)  { fprintf(stderr, " %10s", ftit[q]); }
      }
    fprintf(stderr, "\n");
    
    fprintf(stderr, "%6s", "------");
    for (int32_t m = 0; m < 2; m++)
      { fprintf(stderr, " ");
        for (int32_t q = 0; q < nq; q++)  { fprintf(stderr, " %10s", "----------"); }
      }
    fprintf(stderr, "\n");

    fprintf(stderr, "%6s", "year");
    for (int32_t q = 0; q < nq; q++)  { fprintf(stderr, " %10s", ftit[q]); }
    fprintf(stderr, "\n");
    
    fprintf(stderr, "%6s", "------");
    for (int32_t q = 0; q < nq; q++) { fprintf(stderr, " %10s", "----------"); }
    fprintf(stderr, "\n");
