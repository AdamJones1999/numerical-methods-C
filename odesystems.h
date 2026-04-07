
void orbitalmotion(double drdt_t[], double r_t[], double t, uint Nr);


void orbitalburn(double drdt_t[], double r_t[], double t, uint Nr);


int Schrodinger1DVariableE(double drdx[], double x, double rcol[], uint Nr, void *params);


int Schrodinger1D_boundary_val(double r_last[], double *E, uint Nr);