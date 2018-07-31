/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus, Rajesh Gandham

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

#include "agmg.h"

void gmresUpdate(dlong Nrows,
                 dfloat *x,
                 dfloat **V,
                 dfloat *H,
                 dfloat *s,
                 int Niter,
                 int maxIt){

  dfloat *y = (dfloat *) calloc(Niter, sizeof(dfloat));

  for(int k=Niter-1; k>=0; --k){
    y[k] = s[k];

    for(int m=k+1; m<maxIt; ++m)
      y[k] -= H[k + m*(maxIt+1)]*y[m];

    y[k] /= H[k + k*(maxIt+1)];
  }

  for(int j=0; j<Niter; ++j){
    for(dlong n=0; n<Nrows; ++n)
      x[n] += y[j]*V[j][n];
  }

  free(y);
}

void gmresUpdate(parAlmond_t *parAlmond, dlong Nrows,
                 occa::memory o_x,
                 occa::memory *o_V,
                 dfloat *H,
                 dfloat *s,
                 int Niter,
                 int maxIt){

  dfloat *y = (dfloat *) calloc(Niter+1, sizeof(dfloat));

  for(int k=Niter-1; k>=0; --k){
    y[k] = s[k];

    for(int m=k+1; m<Niter; ++m)
      y[k] -= H[k + m*(maxIt+1)]*y[m];

    y[k] /= H[k + k*(maxIt+1)];
  }

  for(int j=0; j<Niter; ++j){
    vectorAdd(parAlmond, Nrows, y[j], o_V[j], 1.0, o_x);
  }

  free(y);
}

void pgmres(parAlmond_t *parAlmond,
           int maxIt,
           dfloat tol){

  csr *A = parAlmond->levels[0]->A;

  const dlong m = A->Nrows;
  // const dlong n = A->Ncols;

  parAlmond->ktype = GMRES;

  // use parAlmond's buffers
  dfloat *r = parAlmond->levels[0]->rhs;
  dfloat *z = parAlmond->levels[0]->x;

  // initial residual
  dfloat nbLocal = innerProd(m, r, r);
  dfloat nb = 0;
  MPI_Allreduce(&nbLocal,&nb,1,MPI_DFLOAT,MPI_SUM,agmg::comm);
  nb = sqrt(nb);

  //    x = 0;
  dfloat *x  = (dfloat *) calloc(m,sizeof(dfloat));
  setVector(m, x, 0.0);

  //sanity check
  if (nb<=tol) {
    for (dlong i=0;i<m;i++)
      parAlmond->levels[0]->x[i] = x[i];

    free(x); 
    return;
  }

  // M r = b - A*x0
  if(parAlmond->options.compareArgs("PARALMOND CYCLE", "KCYCLE")) {
    kcycle(parAlmond, 0);
  } else if(parAlmond->options.compareArgs("PARALMOND CYCLE", "VCYCLE")) {
    vcycle(parAlmond, 0);
  } else {
    for (dlong k=0;k<m;k++)
      z[k] = r[k];  
  }
  for (dlong k=0;k<m;k++)
    r[k] = z[k];

  dfloat nr = innerProd(m, r, r);
  nr = sqrt(nr);

  dfloat *s = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  s[0] = nr;

  dfloat **V = (dfloat **) calloc(maxIt,sizeof(dfloat *));

  for(int i=0; i<maxIt; ++i){
    V[i] = (dfloat *) calloc(m, sizeof(dfloat)); //TODO this is way too much memory if maxit is large
  }

  // V(:,0) = r/nr
  vectorAdd(m, (1./nr), r,  0., V[0]);

  dfloat *H = (dfloat*) calloc((maxIt+1)*(maxIt+1), sizeof(dfloat));
  dfloat *J = (dfloat*) calloc(4*maxIt, sizeof(dfloat));

  dfloat *Av = (dfloat *) calloc(m, sizeof(dfloat));
  dfloat *w  = (dfloat *) calloc(m, sizeof(dfloat));

  int Niter=0;

  for(int i=0; i<maxIt; i++){

    Niter = i+1;
    // Av = A*V(:.i)
    axpy(A, 1.0, V[i], 0.0, Av,parAlmond->nullSpace,parAlmond->nullSpacePenalty);

    // M w = A vi
    for (dlong k=0;k<m;k++)
      r[k] = Av[k];
    if(parAlmond->options.compareArgs("PARALMOND CYCLE", "KCYCLE")) {
      kcycle(parAlmond, 0);
    } else if(parAlmond->options.compareArgs("PARALMOND CYCLE", "VCYCLE")) {
      vcycle(parAlmond, 0);
    } else {
      for (dlong k=0;k<m;k++)
        z[k] = r[k];  
    }
    for (dlong k=0;k<m;k++)
      w[k] = z[k];

    for(int k=0; k<=i; ++k){
      dfloat hkiLocal = innerProd(m, w, V[k]);
      dfloat hki = 0.;
      MPI_Allreduce(&hkiLocal,&hki,1,MPI_DFLOAT,MPI_SUM,agmg::comm);

      // w = w - hki*V[k]
      vectorAdd(m, -hki, V[k], 1.0, w);

      // H(k,i) = hki
      H[k + i*(maxIt+1)] = hki;
    }

    dfloat wdotwLocal = innerProd(m, w, w);
    dfloat wdotw = 0.;
    MPI_Allreduce(&wdotwLocal,&wdotw,1,MPI_DFLOAT,MPI_SUM,agmg::comm);

    H[i+1 + i*(maxIt+1)] = sqrt(wdotw);

    for(int k=0; k<i; ++k){
      dfloat h1 = H[k +     i*(maxIt+1)];
      dfloat h2 = H[k + 1 + i*(maxIt+1)];

      H[k +     i*(maxIt+1)] = J[4*k    ]*h1 + J[4*k + 2]*h2;
      H[k + 1 + i*(maxIt+1)] = J[4*k + 1]*h1 + J[4*k + 3]*h2;
    }

    dfloat h1 = H[i + i*(maxIt+1)];
    dfloat h2 = H[i + 1 + i*(maxIt+1)];
    dfloat hr = sqrt(h1*h1 + h2*h2);

    H[i   +  i*(maxIt+1)] = hr;
    H[i+1 +  i*(maxIt+1)] = 0.;

    dfloat ct = h1/hr;
    dfloat st = h2/hr;
    J[4*i    ] =  ct;     J[4*i + 2] = st;
    J[4*i + 1] = -st;     J[4*i + 3] = ct;

    dfloat s1 = s[i];
    dfloat s2 = s[i+1];

    s[i  ] =  ct*s1 + st*s2;
    s[i+1] = -st*s1 + ct*s2;

    if(fabs(s[i+1]) < tol) break;

    if(i < maxIt-1){
      dfloat wdotwLocal = innerProd(m, w, w);
      dfloat wdotw = 0.;
      MPI_Allreduce(&wdotwLocal,&wdotw,1,MPI_DFLOAT,MPI_SUM,agmg::comm);

      dfloat nw = sqrt(wdotw);

      // V(:,i+1) = w/nw
      vectorAdd(m,1./nw, w, 0.0, V[i+1]);
    }
  }

  gmresUpdate(m, x, V, H, s, Niter, maxIt);

  //copy result back to parAlmond's x storage
  for (dlong i=0;i<m;i++)
    parAlmond->levels[0]->x[i] = x[i];

  free(x); 
  free(s); free(V);
  free(H); free(J);
  free(Av); free(w);

  if(Niter == maxIt)
    printf("gmres did not converge in given number of iterations\n");
}

//TODO need to link this with MPI
void device_pgmres(parAlmond_t *parAlmond,
           int maxIt,
           dfloat tol){

  hyb* A = parAlmond->levels[0]->deviceA;

  const dlong m = A->Nrows;
  // const dlong n = A->Ncols;

  // use parAlmond's buffers
  occa::memory &o_r = parAlmond->levels[0]->o_rhs;
  occa::memory &o_z = parAlmond->levels[0]->o_x;

  // initial residual
  dfloat nbLocal = innerProd(parAlmond, m, o_r, o_r);
  dfloat nb = 0;
  MPI_Allreduce(&nbLocal,&nb,1,MPI_DFLOAT,MPI_SUM,agmg::comm);
  nb = sqrt(nb);

  dfloat *dummy = (dfloat*) calloc(m, sizeof(dfloat));
  occa::memory  o_x = parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  occa::memory  o_Av= parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  occa::memory  o_w = parAlmond->device.malloc(m*sizeof(dfloat), dummy);

  //sanity check
  if (nb<=tol) {
    parAlmond->levels[0]->o_x.copyFrom(o_x);
    printf("Almond PGMRES iter %d, res = %g\n", 0, nb);
    o_x.free(); o_Av.free(); o_w.free();
    return;
  }

  // M r = b - A*x0
  if(parAlmond->options.compareArgs("PARALMOND CYCLE", "KCYCLE")) {
    device_kcycle(parAlmond, 0);
  } else if(parAlmond->options.compareArgs("PARALMOND CYCLE", "VCYCLE")) {
    device_vcycle(parAlmond, 0);
  } else {
    o_z.copyFrom(o_r);
  }
  o_r.copyFrom(o_z);

  dfloat nrLocal = innerProd(parAlmond, m, o_r, o_r);
  dfloat nr = 0;
  MPI_Allreduce(&nrLocal,&nr,1,MPI_DFLOAT,MPI_SUM,agmg::comm);
  nr = sqrt(nr);

  dfloat *s = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  s[0] = nr;

  occa::memory *o_V = (occa::memory *) calloc(maxIt, sizeof(occa::memory));
  for(int i=0; i<maxIt; ++i){
    o_V[i] = parAlmond->device.malloc(m*sizeof(dfloat), dummy);
  }
  free(dummy);

  // V(:,0) = r/nr
  vectorAdd(parAlmond, m, (1./nr), o_r, 0., o_V[0]);

  dfloat *H = (dfloat *) calloc((maxIt+1)*(maxIt+1), sizeof(dfloat));
  dfloat *J = (dfloat *) calloc(4*maxIt, sizeof(dfloat));

  int Niter = 0;

  int i;
  for(i=0; i<maxIt; i++){

    Niter = i+1;

    // r = A*V(:.i)
    axpy(parAlmond, A, 1.0, o_V[i], 0.0, o_r,parAlmond->nullSpace,parAlmond->nullSpacePenalty);

    // M w = A vi
    if(parAlmond->options.compareArgs("PARALMOND CYCLE", "KCYCLE")) {
      device_kcycle(parAlmond, 0);
    } else if(parAlmond->options.compareArgs("PARALMOND CYCLE", "VCYCLE")) {
      device_vcycle(parAlmond, 0);
    } else {
      o_z.copyFrom(o_r);
    }

    for(int k=0; k<=i; ++k){
      dfloat hkiLocal = innerProd(parAlmond, m, o_z, o_V[k]);
      dfloat hki = 0.;
      MPI_Allreduce(&hkiLocal,&hki,1,MPI_DFLOAT,MPI_SUM,agmg::comm);

      // w = w - hki*V[k]
      vectorAdd(parAlmond, m, -hki, o_V[k], 1.0, o_z);

      // H(k,i) = hki
      H[k + i*(maxIt+1)] = hki;
    }

    dfloat nwLocal = innerProd(parAlmond, m, o_z, o_z);
    dfloat nw = 0.;
    MPI_Allreduce(&nwLocal,&nw,1,MPI_DFLOAT,MPI_SUM,agmg::comm);
    nw = sqrt(nw);
    H[i+1 + i*(maxIt+1)] = nw;

    for(int k=0; k<i; ++k){
      dfloat h1 = H[k +     i*(maxIt+1)];
      dfloat h2 = H[k + 1 + i*(maxIt+1)];

      H[k +     i*(maxIt+1)] = J[4*k    ]*h1 + J[4*k + 2]*h2;
      H[k + 1 + i*(maxIt+1)] = J[4*k + 1]*h1 + J[4*k + 3]*h2;
    }

    dfloat h1 = H[i + i*(maxIt+1)];
    dfloat h2 = H[i + 1 + i*(maxIt+1)];
    dfloat hr = sqrt(h1*h1 + h2*h2);

    H[i   +  i*(maxIt+1)] = hr;
    H[i+1 +  i*(maxIt+1)] = 0;

    dfloat ct = h1/hr;
    dfloat st = h2/hr;
    J[4*i    ] =  ct;     J[4*i + 2] = st;
    J[4*i + 1] = -st;     J[4*i + 3] = ct;

    dfloat s1 = s[i];
    dfloat s2 = s[i+1];

    s[i  ] =  ct*s1 + st*s2;
    s[i+1] = -st*s1 + ct*s2;

    if(fabs(s[i+1]) < tol) break;

    if(i < maxIt-1){
      // V(:,i+1) = w/nw
      vectorAdd(parAlmond, m, 1./nw, o_z, 0.0, o_V[i+1]);
    }
  }

  gmresUpdate(parAlmond, m, o_x, o_V, H, s, Niter, maxIt);

  //copy result back to parAlmond's x storage
  parAlmond->levels[0]->o_x.copyFrom(o_x);

  printf("Almond PGMRES iter %d, res = %g\n", Niter, fabs(s[i+1]));

  if(Niter == maxIt)
    printf("gmres did not converge in given number of iterations \n");

  for(int i=0; i<maxIt; ++i)
    o_V[i].free();
  free((void*)o_V);

  free(s); 
  free(H); free(J);

  o_Av.free();
  o_w.free();
  o_x.free();
}
