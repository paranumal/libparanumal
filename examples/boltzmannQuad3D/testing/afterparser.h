#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <complex.h>

#define dfloat double
#define iint int

typedef struct {
  dfloat *velocity;
  dfloat *vorticity;
  dfloat *density;
  dfloat *q;
  dfloat *jacobian;
}
