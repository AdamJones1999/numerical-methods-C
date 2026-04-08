
void orbitalmotion(double drdt_t[], double r_t[], double t, uint Nr);


void orbitalburn(double drdt_t[], double r_t[], double t, uint Nr);


void Schrodinger1DVariableE(double drdx[], double rcol[], double x, uint Nr, void *E);


void Schrodinger1D_boundary_val(double r_last[], double *E, uint Nr);