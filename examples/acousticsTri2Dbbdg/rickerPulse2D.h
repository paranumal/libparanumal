
//Ricker pulse
#define ricker(t, f) \
 (1-2*OCCA_PI*OCCA_PI*f*f*t*t)*occaExp(-OCCA_PI*OCCA_PI*f*f*t*t) 

//integrated Ricker pulse
#define intRicker(t, f) \
  t*occaExp(-OCCA_PI*OCCA_PI*f*f*t*t) 

#define acousticsPointSource2D(x, y, t, f, c, u, v, p)  \
{ \
  dfloat r = sqrt(x*x+y*y); \
                              \
  p = ricker((t-r/c),f)/(4*OCCA_PI*c*c*r); \
  u = x*(intRicker((t-r/c),f)/r + ricker((t-r/c),f)/c)/(4*OCCA_PI*c*c*r*r); \
  v = y*(intRicker((t-r/c),f)/r + ricker((t-r/c),f)/c)/(4*OCCA_PI*c*c*r*r); \
}