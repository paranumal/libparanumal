
//Ricker pulse
#define ricker(t, f) \
 (1-2*OCCA_PI*OCCA_PI*f*f*t*t)*occaExp(-OCCA_PI*OCCA_PI*f*f*t*t)

//integrated Ricker pulse
#define intRicker(t, f) \
  t*occaExp(-OCCA_PI*OCCA_PI*f*f*t*t)

#define acousticsPointSource3D(x, y, z, t, f, c, u, v, w, p)  \
{ \
  dfloat r = sqrt(x*x+y*y+z*z); \
                              \
  p = ricker((t-r/c),f)/(4*OCCA_PI*c*c*r); \
  u = x*(intRicker((t-r/c),f)/r + ricker((t-r/c),f)/c)/(4*OCCA_PI*c*c*r*r); \
  v = y*(intRicker((t-r/c),f)/r + ricker((t-r/c),f)/c)/(4*OCCA_PI*c*c*r*r); \
  w = z*(intRicker((t-r/c),f)/r + ricker((t-r/c),f)/c)/(4*OCCA_PI*c*c*r*r); \
}