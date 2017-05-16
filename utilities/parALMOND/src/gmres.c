#include "parAlmond.h"

void gmresUpdate(iint Nrows,
                 dfloat *x,
            	   dfloat **V,
            	   dfloat *H,
            	   dfloat *s,
            	   iint end,
            	   iint maxIt){

  dfloat *y = (dfloat *) calloc(end, sizeof(dfloat));

  for(iint k=end-1; k>=0; --k){
    y[k] = s[k];

    for(iint m=k+1; m<maxIt; ++m)
      y[k] -= H[k + m*(maxIt+1)]*y[m];

    y[k] /= H[k + k*(maxIt+1)];
  }

  for(iint j=0; j<end; ++j){
    for(iint n=0; n<Nrows; ++n)
      x[n] += y[j]*V[j][n];
  }
}

void gmresUpdate(parAlmond_t *parAlmond, iint Nrows,
     occa::memory o_x,
	   occa::memory *o_V,
	   dfloat *H,
	   dfloat *s,
	   iint end,
	   iint maxIt){

  dfloat *y = (dfloat *) calloc(end+1, sizeof(dfloat));

  for(iint k=end-1; k>=0; --k){
    y[k] = s[k];

    for(iint m=k+1; m<end; ++m)
  	  y[k] -= H[k + m*(maxIt+1)]*y[m];

    y[k] /= H[k + k*(maxIt+1)];
  }

  for(iint j=0; j<end; ++j){
    vectorAdd(parAlmond, Nrows, y[j], o_V[j], 1.0, o_x);
  }
}

//TODO need to link this with MPI
void gmres(parAlmond_t *parAlmond,
     csr *A,
     dfloat *b,
     dfloat *x,
     iint maxIt,
     dfloat tol){

  iint m = A->Nrows;
  iint n = A->Ncols;

  parAlmond->ktype = GMRES;

  // initial residual
  dfloat nb = norm(m, b);

  dfloat *r = (dfloat *) calloc(m,sizeof(dfloat));

  // M r = b - A*x0
  //solve(parAlmond, b, r);
  for (iint k=0;k<m;k++)
    r[k] = b[k];

  dfloat nr = norm(m, r);

  dfloat *s = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  s[0] = nr;

  dfloat **V = (dfloat **) calloc(maxIt,sizeof(dfloat *));

  for(iint i=0; i<maxIt; ++i){
    V[i] = (dfloat *) calloc(m, sizeof(dfloat)); //TODO this is way too much memory if maxit is large
  }

  // V(:,0) = r/nr
  vectorAdd(m, (1./nr), r,  0., V[0]);

  dfloat *H = (dfloat*) calloc((maxIt+1)*(maxIt+1), sizeof(dfloat));
  dfloat *J = (dfloat*) calloc(4*maxIt, sizeof(dfloat));

  dfloat *resVec = (dfloat *) calloc(maxIt, sizeof(dfloat));

  resVec[0] = nr;
  dfloat *Av = (dfloat *) calloc(m, sizeof(dfloat));
  dfloat *w  = (dfloat *) calloc(m, sizeof(dfloat));

  iint end;

  for(iint i=0; i<maxIt; i++){

    end = i+1;
    // Av = A*V(:.i)
    axpy(A, 1.0, V[i], 0.0, Av);

    // M w = A vi
    //solve(parAlmond, Av, w);
    for (iint k=0;k<m;k++)
      w[k] = Av[k];

    for(iint k=0; k<=i; ++k){
      dfloat hki = innerProd(m, w, V[k]);

      // w = w - hki*V[k]
      vectorAdd(m, -hki, V[k], 1.0, w);

      // H(k,i) = hki
      H[k + i*(maxIt+1)] = hki;
    }

    H[i+1 + i*(maxIt+1)] = norm(m,w);


    for(iint k=0; k<i; ++k){
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

    resVec[i+1] = fabs(s[i+1]);

    if(fabs(s[i+1]) < tol) break;
    
    if(i < maxIt-1){
      dfloat nw = norm(m,w);

      // V(:,i+1) = w/nw
      vectorAdd(m,1./nw, w, 0.0, V[i+1]);
    }
  }

  gmresUpdate(m, x, V, H, s, end, maxIt);

  if(end == maxIt)
    printf("gmres did not converge in given number of iterations\n");
}

//TODO need to link this with MPI
void gmres(parAlmond_t *parAlmond,
     hyb *A,
     occa::memory o_b,
     occa::memory o_x,
     iint maxIt,
     dfloat tol){


  iint m = A->Nrows;
  iint n = A->Ncols;

  iint sz = m*sizeof(dfloat);

  parAlmond->ktype = GMRES;

  // initial residual
  dfloat nb = innerProd(parAlmond, m, o_b, o_b);
  nb = sqrt(nb);

  dfloat *dummy = (dfloat*) calloc(m, sizeof(dfloat));

  occa::memory  o_r = parAlmond->device.malloc(sz, dummy);
  occa::memory o_Av = parAlmond->device.malloc(sz, dummy);
  occa::memory  o_w = parAlmond->device.malloc(sz, dummy);

  // M r = b - A*x0
  //solve(parAlmond, o_b, o_r);
  o_r.copyFrom(o_b);


  dfloat nr = innerProd(parAlmond, m, o_r, o_r);
  nr = sqrt(nr);

  dfloat *s = (dfloat *) calloc(maxIt+1, sizeof(dfloat));
  s[0] = nr;

  occa::memory *o_V = (occa::memory *) calloc(maxIt, sizeof(occa::memory));

  for(iint i=0; i<maxIt; ++i){
    o_V[i] = parAlmond->device.malloc(sz, dummy);
  }

  // V(:,0) = r/nr
  vectorAdd(parAlmond, m, (1./nr), o_r, 0., o_V[0]);

  dfloat *H = (dfloat *) calloc((maxIt+1)*(maxIt+1), sizeof(dfloat));
  dfloat *J = (dfloat *) calloc(4*maxIt, sizeof(dfloat));
  dfloat *resVec = (dfloat *) calloc(maxIt, sizeof(dfloat));

  resVec[0] = nr;

  iint end = 0;

  for(iint i=0; i<maxIt; i++){

    end = i+1;

    // Av = A*V(:.i)
    axpy(parAlmond, A, 1.0, o_V[i], 0.0, o_Av);

    // M w = A vi
    //solve(parAlmond, o_Av, o_w);
    o_w.copyFrom(o_Av);

    for(iint k=0; k<=i; ++k){
      dfloat hki = innerProd(parAlmond, m, o_w, o_V[k]);

      // w = w - hki*V[k]
      vectorAdd(parAlmond, m, -hki, o_V[k], 1.0, o_w);

      // H(k,i) = hki
      H[k + i*(maxIt+1)] = hki;
    }

    dfloat nw = innerProd(parAlmond, m, o_w, o_w);
    nw = sqrt(nw);
    H[i+1 + i*(maxIt+1)] = nw;


    for(iint k=0; k<i; ++k){
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

    resVec[i+1] = fabs(s[i+1]);

    if(fabs(s[i+1]) < tol) break;

    if(i < maxIt-1){
      // V(:,i+1) = w/nw
      vectorAdd(parAlmond, m, 1./nw, o_w, 0.0, o_V[i+1]);
    }
  }

  gmresUpdate(parAlmond, m, o_x, o_V, H, s, end, maxIt);

  if(end == maxIt)
    printf("gmres did not converge in given number of iterations \n");
  
  for(iint i=0; i<maxIt; ++i)
    o_V[i].free();

  o_Av.free();
  o_w.free();
  o_r.free();
}

