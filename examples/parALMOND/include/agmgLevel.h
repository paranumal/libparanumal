
void agmgAx        (void **args, dfloat *x, dfloat *Ax);
void agmgCoarsen   (void **args, dfloat *r, dfloat *Rr);
void agmgProlongate(void **args, dfloat *x, dfloat *Px);
void agmgSmooth    (void **args, dfloat *rhs, dfloat *x, bool x_is_zero);

void device_agmgAx        (void **args, occa::memory &o_x, occa::memory &o_Ax);
void device_agmgCoarsen   (void **args, occa::memory &o_r, occa::memory &o_Rr);
void device_agmgProlongate(void **args, occa::memory &o_x, occa::memory &o_Px);
void device_agmgSmooth    (void **args, occa::memory &o_r, occa::memory &o_x, bool x_is_zero);

void setupSmoother(agmgLevel *level, SmoothType s);
void setupExactSolve(parAlmond_t *parAlmond, agmgLevel *level);
