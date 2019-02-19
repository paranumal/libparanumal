#include "adaptive.h"

void get_operators(int N, occa::device &device, occa::memory &o_r,
                   occa::memory &o_w, occa::memory &o_D, occa::memory &o_Ib,
                   occa::memory &o_It, occa::memory &o_Pb, occa::memory &o_Pt)
{
  const int Nq = N + 1;

  long double *lr = (long double *)asd_malloc_aligned(Nq * sizeof(long double));
  long double *lw = (long double *)asd_malloc_aligned(Nq * sizeof(long double));
  long double *lV = (long double *)asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lD = (long double *)asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lM = (long double *)asd_malloc_aligned(Nq * Nq * sizeof(long double));

  long double *lrb = (long double *)asd_malloc_aligned(Nq * sizeof(long double));
  long double *lrt = (long double *)asd_malloc_aligned(Nq * sizeof(long double));

  long double *lIb = (long double *)asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lIt = (long double *)asd_malloc_aligned(Nq * Nq * sizeof(long double));

  long double *lPb = (long double *)asd_malloc_aligned(Nq * Nq * sizeof(long double));
  long double *lPt = (long double *)asd_malloc_aligned(Nq * Nq * sizeof(long double));

  dfloat_t *I = (dfloat_t*)asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));
  for (int n = 0; n < Nq * Nq; ++n)
    I[n] = 0;
  for (int n = 0; n < Nq; ++n)
    I[n + n * Nq] = 1;

  asd_jacobi_gauss_lobatto_quadrature(0, 0, N, lr, lw);
  asd_jacobi_p_vandermonde(0, 0, N, Nq, lr, lV);
  asd_jacobi_p_differentiation(0, 0, N, Nq, lr, lV, lD);
  asd_jacobi_p_mass(N, lV, lM);

  for (int n = 0; n < Nq; ++n)
  {
    lrb[n] = (lr[n] - 1) / 2;
    lrt[n] = (lr[n] + 1) / 2;
  }

  asd_jacobi_p_interpolation(0, 0, N, Nq, lrb, lV, lIb);
  asd_jacobi_p_interpolation(0, 0, N, Nq, lrt, lV, lIt);

  asd_jacobi_p_h_project(N, 0.5L, lV, lIb, lM, lPb);
  asd_jacobi_p_h_project(N, 0.5L, lV, lIt, lM, lPt);

  dfloat_t *r = (dfloat_t*)asd_malloc_aligned(Nq * sizeof(dfloat_t));
  dfloat_t *w = (dfloat_t*)asd_malloc_aligned(Nq * sizeof(dfloat_t));

  for (int n = 0; n < Nq; ++n)
  {
    r[n] = (dfloat_t)lr[n];
    w[n] = (dfloat_t)lw[n];
  }

  o_r = device.malloc(Nq * sizeof(dfloat_t), r);
  o_w = device.malloc(Nq * sizeof(dfloat_t), w);

  dfloat_t *DT = (dfloat_t*)asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));

  for (int i = 0; i < Nq; ++i)
    for (int j = 0; j < Nq; ++j)
      DT[j * Nq + i] = (dfloat_t)lD[i * Nq + j];

  o_D = device.malloc(Nq * Nq * sizeof(dfloat_t), DT);

  dfloat_t *Ib = (dfloat_t*)asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));
  dfloat_t *It = (dfloat_t*)asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));

  for (int i = 0; i < Nq; ++i)
  {
    for (int j = 0; j < Nq; ++j)
    {
      Ib[j * Nq + i] = (dfloat_t)lIb[i * Nq + j];
      It[j * Nq + i] = (dfloat_t)lIt[i * Nq + j];
    }
  }

  o_Ib = device.malloc(Nq * Nq * sizeof(dfloat_t), Ib);
  o_It = device.malloc(Nq * Nq * sizeof(dfloat_t), It);

  dfloat_t *Pb = (dfloat_t*)asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));
  dfloat_t *Pt = (dfloat_t*)asd_malloc_aligned(Nq * Nq * sizeof(dfloat_t));

  for (int i = 0; i < Nq; ++i)
  {
    for (int j = 0; j < Nq; ++j)
    {
      Pb[j * Nq + i] = (dfloat_t)lPb[i * Nq + j];
      Pt[j * Nq + i] = (dfloat_t)lPt[i * Nq + j];
    }
  }

  printf("Ib = [\n");
  for (int i = 0; i < Nq; ++i){
    for (int j = 0; j < Nq; ++j){
      printf("%17.15lf ", Ib[j*Nq+i]);
    }
    printf("\n");
  }
  printf("];\n");

  printf("It = [\n");
  for (int i = 0; i < Nq; ++i){
    for (int j = 0; j < Nq; ++j){
      printf("%17.15lf ", It[j*Nq+i]);
    }
    printf("\n");
  }
  printf("];\n");

  o_Pb = device.malloc(Nq * Nq * sizeof(dfloat_t), Pb);
  o_Pt = device.malloc(Nq * Nq * sizeof(dfloat_t), Pt);
  
  asd_free_aligned(lr);
  asd_free_aligned(lw);
  asd_free_aligned(lV);
  asd_free_aligned(lD);
  asd_free_aligned(lM);
  asd_free_aligned(lrb);
  asd_free_aligned(lrt);
  asd_free_aligned(lIb);
  asd_free_aligned(lIt);
  asd_free_aligned(lPb);
  asd_free_aligned(lPt);

  asd_free_aligned(r);
  asd_free_aligned(w);
  asd_free_aligned(I);
  asd_free_aligned(DT);
  asd_free_aligned(Ib);
  asd_free_aligned(It);
  asd_free_aligned(Pb);
  asd_free_aligned(Pt);
}
