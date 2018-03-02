#include <stdio.h>
#include <math.h>
#include <complex.h>

#define dfloat double
#define iint int

#define BUFSIZ 1000

typedef struct {
  dfloat *velocity;
  dfloat *vorticity;
  dfloat *density;
  dfloat *q;
  dfloat *jacobian;
} data;

void parse_vtu(char *filename, data *parsed);
