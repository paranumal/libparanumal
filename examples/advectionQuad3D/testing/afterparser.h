#ifndef afterparser
#define afterparser

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

#define dfloat double
#define iint int

typedef struct {
  iint points;
  iint elements;
  dfloat *x;
  dfloat *y;
  dfloat *z;
  iint *mrab_levels;
  dfloat *vel_x;
  dfloat *vel_y;
  dfloat *vel_z;
  dfloat *vort_x;
  dfloat *vort_y;
  dfloat *vort_z;
  dfloat *density;
  dfloat *q;
  dfloat *jacobian;
  dfloat *rad_vort;
  iint *connectivity; //too lazy to generate this myself...
  
  //error terms.  Should make a new struct when there are more
  dfloat *mismatch;
  dfloat *difference;
} data;

void parse_vtu(char *filename, data *parsed);

void send_vtu(char *filename, data *parsed);

dfloat compute_l2(data *ref, data *test);

#endif
