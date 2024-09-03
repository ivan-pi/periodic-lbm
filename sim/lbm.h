
// These are the callback types

// The name of the procedures will be the name of the module.
// For instance if the module is libslbm.so, the procedures
// should be named
//
// f_slbminit()
// f_slbmstep()
// f_slbminit()
// f_slbminit()

void *siminit(
	int nx, int ny, double dt,
	double *rho, double *u, double *sigma, void *params);

void simstep(void *sim, double omega);
void simvars(void *sim, double *rho, double *u);
void simfree(void *sim);




