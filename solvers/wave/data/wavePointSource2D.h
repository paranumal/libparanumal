
/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#define ellipticForcing2D(x, y, lambda, f) \
  { \
  f = 0.;                                     \
  }


/* Homogeneous Dirichlet boundary condition   */
#define ellipticDirichletCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = 0.f;   \
    uxB = uxM;   \
    uyB = uyM;   \
  }

/* Homogeneous Neumann boundary condition   */
#define ellipticNeumannCondition2D(x,y,nx,ny,uM,uxM,uyM,uB,uxB,uyB)  \
  {              \
    uB  = uM;    \
    uxB = 0.f;   \
    uyB = 0.f;   \
  }

#define waveForcingFunction2D(x, y, sigma, omega, f) \
  {                                            \
    /*    f = exp(-80.*(x*x+(y-.5)*(y-.5)));    */  \
    f = 0.;                                         \
  }

#define waveRicker2D(f, t, ricker, intRicker, drickerdt, drickerdt2)    \
  {                                                                     \
    dfloat fac = M_PI*M_PI*f*f;                                         \
    dfloat sc  = exp(-fac*t*t);                                          \
    intRicker  = t*sc;                                                  \
    ricker     = (1.f-2.f*fac*t*t)*sc;                                  \
    drickerdt  =  2.f*fac*t*(2.f*fac*t*t - 3.f)*sc;                     \
    drickerdt2 = -2.f*fac*(4.f*fac*fac*t*t*t*t-12.f*fac*t*t+3.f)*sc;    \
  }


#define waveInitialConditionsFunction2D(t, x, y, d, p)  \
  {                                                 \
    d = 0.;                                         \
    p = 0.;                                         \
  }

//    d    =  intRicker/(4.f*M_PI*R);                                   
//    dddx = -dRdx*ricker/(4.f*M_PI*R) - intRicker*dRdx/(4.f*M_PI*R*R);   
//    dddy = -dRdy*ricker/(4.f*M_PI*R) - intRicker*dRdy/(4.f*M_PI*R*R);   

#define waveSurfaceSource2D(patch, tin, xsource, ysource, fsource, xn, yn, p, dpdx, dpdy, d, dddx, dddy) \
  {                                                                     \
    dfloat ricker, intRicker, drickerdt, drickerdt2;                    \
    dfloat xs = xn-xsource;                                             \
    dfloat ys = yn-ysource;                                             \
    dfloat R = sqrt(xs*xs+ys*ys);                                       \
    dfloat dRdx = xs/R;                                                 \
    dfloat dRdy = ys/R;                                                 \
    dfloat tMR = tin-R-.5f;                                             \
                                                                        \
    waveRicker2D(fsource, (tMR), ricker, intRicker, drickerdt, drickerdt2); \
                                                                        \
    p    = ricker/(4.*M_PI*R);                                          \
    dpdx = -dRdx*drickerdt/(4.*M_PI*R) - ricker*dRdx/(4.*M_PI*R*R);     \
    dpdy = -dRdy*drickerdt/(4.*M_PI*R) - ricker*dRdy/(4.*M_PI*R*R);     \
                                                                        \
    /* d = -div U = dpdt*/                                              \
    /* CHANGED THIS TO integral rather than derivative */               \
                                                                        \
    d    = drickerdt/(4.f*M_PI*R);                                      \
    dddx = -dRdx*drickerdt2/(4.f*M_PI*R) - drickerdt*dRdx/(4.f*M_PI*R*R); \
    dddy = -dRdy*drickerdt2/(4.f*M_PI*R) - drickerdt*dRdy/(4.f*M_PI*R*R); \
                                                                        \
    if(patch==-1){                                                       \
      /* not 100 percent sure */                                        \
      p *= -1.f; dpdx *= -1.f; dpdy *= -1.f;                            \
      d *= -1.f; dddx *= -1.f; dddy *= -1.f;                            \
    }                                                                   \
    if(patch==2){                                                       \
      /* not 100 percent sure */                                        \
      /* dirichlet p+ = -p-, dpdx+ = dpdx- */                           \
      p *= -2.f; dpdx *= -2.f; dpdy *= -2.f;                            \
      d *= -2.f; dddx *= -2.f; dddy *= -2.f;                            \
    }                                                                   \
    /*    printf("tX=%e,%e, %e == p=%e,%e,%e === d=%e,%e,%e\n", tin, xs, ys, p, dpdx, dpdy, d, dddx, dddy); */ \
  }
    
