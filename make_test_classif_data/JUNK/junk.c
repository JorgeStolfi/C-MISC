/* Last edited on 2010-05-22 06:55:09 by stolfi */

/* PRBOBLEM SETUP PROCS */
  
void problem_def_saturn(int NC, problem_t *P)
  { if (NC == 0) { NC = 2; }
    demand(NC == 2, "invalid num of classes for this problem");
    P->NC = NC;
    P->NA = 2;
    P->sampler = &problem_smp_saturn;
  }

void problem_def_petals(int NC, problem_t *P)
  { if (NC == 0) { NC = 4; }
    demand(NC == 4, "invalid num of classes for this problem");
    P->NC = NC;
    P->NA = 2;
    P->sampler = &problem_smp_petals;
    
  }

void problem_def_vessel(int NC, problem_t *P)
  { if (NC == 0) { NC = 3; }
    demand(NC == 3, "invalid num of classes for this problem");
    P->NC = NC;
    P->NA = 2;
    P->sampler = &problem_smp_vessel;
  }

void problem_def_mballs(int NC, problem_t *P)
  { if (NC == 0) { NC = 2; }
    demand(NC >= 3, "invalid num of classes for this problem");
    P->NC = NC;
    P->NA = 2;
    P->sampler = &problem_smp_mballs;
  }

void problem_def_shells(int NC, problem_t *P)
  { if (NC == 0) { NC = 2; }
    demand(NC >= 2, "invalid num of classes for this problem");
    P->NC = NC;
    P->NA = 2;
    P->sampler = &problem_smp_shells;
  }
