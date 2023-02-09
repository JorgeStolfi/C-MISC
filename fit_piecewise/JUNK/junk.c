
void fpc_build_model(int NP, double Z[], double W[], int NX, double X[], double C[]);
  /* Takes a list {Z[0..NP-1]} of data samples with weights {W[0..NP-1]}, and a linearized {NP} by
    {NX} array {X[0..NP*NX-1]} of the corresponding values of the
    independent variables. Builds a linear model that gives an
    approximation {Y[i]} to each {Z[i]} by a combination of the
    independent variables {X[i,0..NX-1]}.
    
    The model is {Y = X*C} where {Y[0..NP-1]} is the vector of predicted values,
    and {C[0..NX-1]} is a columnn vector of coefficients that are
    determined by the procedure. */
  
void fpc_apply_model(int NP, int NX, double X[], double C[], double Y[]);
  /* Takes a series {Z[0..NP-1]} of data samples, and the coefficients
    {C[0..NX-1]} of a linear predictor. Computes the approximation values
    {Y[0..NP-1]}. See {fpc_build_model} for details. */
   
void fpc_write_model(char *fname, int NX, double C[], double avg, double dev);
  /* Writes the fitted model to file "{fname}", in a machine-readable format. */
  
void fpc_print_model(FILE *wr, int NX, double C[], double avg, double dev);
  /* Prints the fitted model to {wr}. */
