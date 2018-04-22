#include "advectionQuad3D.h"

void rk4_coeffs(solver_t *solver) {
  int Nrk = 5;
  
  solver->rka = (dfloat *) calloc(Nrk,sizeof(dfloat));
  solver->rkb = (dfloat *) calloc(Nrk,sizeof(dfloat));
  solver->rkc = (dfloat *) calloc(Nrk+1,sizeof(dfloat));
  
  dfloat rka[5] = {0.0,
		   -567301805773.0/1357537059087.0 ,
		   -2404267990393.0/2016746695238.0 ,
		   -3550918686646.0/2091501179385.0  ,
		   -1275806237668.0/842570457699.0};
  dfloat rkb[5] = { 1432997174477.0/9575080441755.0 ,
		    5161836677717.0/13612068292357.0 ,
		    1720146321549.0/2090206949498.0  ,
		    3134564353537.0/4481467310338.0  ,
		    2277821191437.0/14882151754819.0};
  dfloat rkc[6] = {0.0  ,
		   1432997174477.0/9575080441755.0 ,
		   2526269341429.0/6820363962896.0 ,
		   2006345519317.0/3224310063776.0 ,
		   2802321613138.0/2924317926251.0,
		   1.};
  solver->Nrk = Nrk;
  memcpy(solver->rka, rka, Nrk*sizeof(dfloat));
  memcpy(solver->rkb, rkb, Nrk*sizeof(dfloat));
  memcpy(solver->rkc, rkc, (Nrk+1)*sizeof(dfloat));
}

void mrab3_coeffs(solver_t *solver) {

  mesh_t *mesh = solver->mesh;
  
  iint Nlevels = mesh->MRABNlevels;

  const iint Nr = 32;
  dfloat complex R[Nr];

  for(iint ind =1; ind <= Nr; ++ind){
    const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
    R[ind-1] = cexp(I*M_PI* theta);
  }

  solver->MRSAAB_A = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  solver->MRSAAB_B = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  solver->MRSAAB_C = (dfloat *) calloc(    Nlevels,sizeof(dfloat));
  solver->MRAB_A   = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  solver->MRAB_B   = (dfloat *) calloc(3*3*Nlevels,sizeof(dfloat));
  solver->MRAB_C   = (dfloat *) calloc(    Nlevels,sizeof(dfloat));

  iint MRABorder = solver->Nrhs;

  for(iint l = 0; l<Nlevels; ++l){
    // MRSAAB coefficients
    dfloat alpha = -solver->tauInv*solver->dt*pow(2,l);
    dfloat h  = solver->dt * pow(2,l);
    //
    for (iint order=0; order<3; ++order){
      // computation of coefficients based on magnitude
      const iint id = order*Nlevels*3 + l*3;
      if(order==0){

	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 +=  h*(cexp(lr) - 1.)/lr;
	  b1 +=  h*(cexp(lr/2.) - 1.)/lr;
	}
	// Full dt coeeficients
	solver->MRSAAB_A[id + 0] = creal(a1)/Nr;
	solver->MRSAAB_A[id + 1] = 0.f;
	solver->MRSAAB_A[id + 2] = 0.f;
	// Half coefficients
	solver->MRSAAB_B[id + 0] = creal(b1)/Nr;
	solver->MRSAAB_B[id + 1] = 0.f;
	solver->MRSAAB_B[id + 2] = 0.f;

	// MRAB coefficients
	solver->MRAB_A[id + 0]   =  h ;
	solver->MRAB_A[id + 1]   =  0.f ;
	solver->MRAB_A[id + 2]   =  0.f ;

	solver->MRAB_B[id+0]     =  h/2. ;
	solver->MRAB_B[id+1]     =  0.f ;
	solver->MRAB_B[id+2]     =  0.f ;
      }

      else if(order==1){

	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;
	double complex a2 = 0. + 0.* I;
	double complex b2 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 +=  h*(-2.*lr + (1.+lr)*cexp(lr) - 1.)/cpow(lr,2);
	  a2 +=  h*(lr - cexp(lr) + 1.)/cpow(lr,2);
	  b1 +=  h*(-1.5*lr + (1.+lr)*cexp(lr/2.) - 1.)/cpow(lr,2);
	  b2 +=  h*(0.5*lr - cexp(lr/2.) + 1.)/cpow(lr,2);
	}
	// Full dt coeeficients
	solver->MRSAAB_A[id + 0] = creal(a1)/Nr;
	solver->MRSAAB_A[id + 1] = creal(a2)/Nr;
	solver->MRSAAB_A[id + 2] = 0.f;
	// Half coefficients
	solver->MRSAAB_B[id + 0] = creal(b1)/Nr;
	solver->MRSAAB_B[id + 1] = creal(b2)/Nr;
	solver->MRSAAB_B[id + 2] = 0.f;


	// MRAB coefficients
	solver->MRAB_A[id + 0]   =  3.*h/2. ;
	solver->MRAB_A[id + 1]   = -1.*h/2. ;
	solver->MRAB_A[id + 2]   =  0.f ;

	solver->MRAB_B[id + 0]   =  5.*h/8. ;
	solver->MRAB_B[id + 1]   = -1.*h/8. ;
	solver->MRAB_B[id + 2]   =   0.f ;
      }

      else{
	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;
	double complex a2 = 0. + 0.* I;
	double complex b2 = 0. + 0.* I;
	double complex a3 = 0. + 0.* I;
	double complex b3 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 += h*(-2.5*lr - 3.*cpow(lr,2) + (1.+cpow(lr,2)+1.5*lr)*cexp(lr) - 1.)/cpow(lr,3);
	  a2 += h*(4.*lr + 3.*cpow(lr,2)- (2.*lr + 2.0)*cexp(lr) + 2.)/cpow(lr,3);
	  a3 +=-h*(1.5*lr + cpow(lr,2)- (0.5*lr + 1.)*cexp(lr) + 1.)/cpow(lr,3);
	  b1 += h*(cexp(lr/2.)- 2.*lr - (15.*cpow(lr,2))/8.f + cpow(lr,2)*cexp(lr/2.) + 3.*lr*cexp(lr/2.)/2. - 1.)/cpow(lr,3);
	  b2 += h*(3.*lr - 2.*cexp(lr/2.0) + 1.25*cpow(lr,2) - 2.*lr*cexp(lr/2.) + 2.)/cpow(lr,3);
	  b3 +=-h*(lr - cexp(lr/2.) + 0.375*cpow(lr,2) - 0.5*lr*cexp(lr/2.) + 1.)/cpow(lr,3);
	}


	// Full dt coeeficients
	solver->MRSAAB_A[id+0] = creal(a1)/Nr;
	solver->MRSAAB_A[id+1] = creal(a2)/Nr;
	solver->MRSAAB_A[id+2] = creal(a3)/Nr;
	// Half coefficients
	solver->MRSAAB_B[id+0] = creal(b1)/Nr;
	solver->MRSAAB_B[id+1] = creal(b2)/Nr;
	solver->MRSAAB_B[id+2] = creal(b3)/Nr;

	// MRAB coefficients
	solver->MRAB_A[id+0]   =  23.*h/12. ;
	solver->MRAB_A[id+1]   = -16.*h/12. ;
	solver->MRAB_A[id+2]   =  5. *h/12. ;

	solver->MRAB_B[id+0]   =  17.*h/24. ;
	solver->MRAB_B[id+1]   = - 7.*h/24. ;
	solver->MRAB_B[id+2]   =   2.*h/24. ;


      }
    }

    // Exponential part
    solver->MRSAAB_C[l]    = exp(alpha);
    solver->MRAB_C[l]      =   h ;
  }
}

void mrab4_coeffs(solver_t *solver) {

  mesh_t *mesh = solver->mesh;
  
  iint Nlevels = mesh->MRABNlevels;

  const iint Nr = 32;
  dfloat complex R[Nr];

  for(iint ind =1; ind <= Nr; ++ind){
    const dfloat theta = (dfloat) (ind - 0.5) / (dfloat) Nr;
    R[ind-1] = cexp(I*M_PI* theta);
  }

  solver->MRSAAB_A = (dfloat *) calloc(4*4*Nlevels,sizeof(dfloat));
  solver->MRSAAB_B = (dfloat *) calloc(4*4*Nlevels,sizeof(dfloat));
  solver->MRSAAB_C = (dfloat *) calloc(    Nlevels,sizeof(dfloat));
  solver->MRAB_A   = (dfloat *) calloc(4*4*Nlevels,sizeof(dfloat));
  solver->MRAB_B   = (dfloat *) calloc(4*4*Nlevels,sizeof(dfloat));
  solver->MRAB_C   = (dfloat *) calloc(    Nlevels,sizeof(dfloat));

  iint MRABorder = solver->Nrhs;

  for(iint l = 0; l<Nlevels; ++l){
    // MRSAAB coefficients
    dfloat alpha = -solver->tauInv*solver->dt*pow(2,l);
    dfloat h  = solver->dt * pow(2,l);
    //
    for (iint order=0; order<4; ++order){
      // computation of coefficients based on magnitude
      const iint id = order*Nlevels*4 + l*4;
      if(order==0){

	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 +=  h*(cexp(lr) - 1.)/lr;
	  b1 +=  h*(cexp(lr/2.) - 1.)/lr;
	}
	// Full dt coeeficients
	solver->MRSAAB_A[id + 0] = creal(a1)/Nr;
	solver->MRSAAB_A[id + 1] = 0.f;
	solver->MRSAAB_A[id + 2] = 0.f;
	solver->MRSAAB_A[id + 3] = 0.f;
	// Half coefficients
	solver->MRSAAB_B[id + 0] = creal(b1)/Nr;
	solver->MRSAAB_B[id + 1] = 0.f;
	solver->MRSAAB_B[id + 2] = 0.f;
	solver->MRSAAB_B[id + 3] = 0.f;

	// MRAB coefficients
	solver->MRAB_A[id + 0]   =  h ;
	solver->MRAB_A[id + 1]   =  0.f ;
	solver->MRAB_A[id + 2]   =  0.f ;
	solver->MRAB_A[id + 3]   =  0.f;

	solver->MRAB_B[id+0]     =  h/2. ;
	solver->MRAB_B[id+1]     =  0.f ;
	solver->MRAB_B[id+2]     =  0.f ;
	solver->MRAB_B[id+3]     =  0.f ;
      }

      else if(order==1){

	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;
	double complex a2 = 0. + 0.* I;
	double complex b2 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 +=  h*(-2.*lr + (1.+lr)*cexp(lr) - 1.)/cpow(lr,2);
	  a2 +=  h*(lr - cexp(lr) + 1.)/cpow(lr,2);
	  b1 +=  h*(-1.5*lr + (1.+lr)*cexp(lr/2.) - 1.)/cpow(lr,2);
	  b2 +=  h*(0.5*lr - cexp(lr/2.) + 1.)/cpow(lr,2);
	}
	// Full dt coeeficients
	solver->MRSAAB_A[id + 0] = creal(a1)/Nr;
	solver->MRSAAB_A[id + 1] = creal(a2)/Nr;
	solver->MRSAAB_A[id + 2] = 0.f;
	solver->MRSAAB_A[id + 3] = 0.f;
	// Half coefficients
	solver->MRSAAB_B[id + 0] = creal(b1)/Nr;
	solver->MRSAAB_B[id + 1] = creal(b2)/Nr;
	solver->MRSAAB_B[id + 2] = 0.f;
	solver->MRSAAB_B[id + 3] = 0.f;

	// MRAB coefficients
	solver->MRAB_A[id + 0]   =  3.*h/2. ;
	solver->MRAB_A[id + 1]   = -1.*h/2. ;
	solver->MRAB_A[id + 2]   =  0.f ;
	solver->MRAB_B[id + 3]   =  0.f ;

	solver->MRAB_B[id + 0]   =  5.*h/8. ;
	solver->MRAB_B[id + 1]   = -1.*h/8. ;
	solver->MRAB_B[id + 2]   =  0.f ;
	solver->MRAB_B[id + 3]   =  0.f ;
      }
      else if (order == 2) {
	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;
	double complex a2 = 0. + 0.* I;
	double complex b2 = 0. + 0.* I;
	double complex a3 = 0. + 0.* I;
	double complex b3 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 += h*(-2.5*lr - 3.*cpow(lr,2) + (1.+cpow(lr,2)+1.5*lr)*cexp(lr) - 1.)/cpow(lr,3);
	  a2 += h*(4.*lr + 3.*cpow(lr,2)- (2.*lr + 2.0)*cexp(lr) + 2.)/cpow(lr,3);
	  a3 +=-h*(1.5*lr + cpow(lr,2)- (0.5*lr + 1.)*cexp(lr) + 1.)/cpow(lr,3);
	  b1 += h*(cexp(lr/2.)- 2.*lr - (15.*cpow(lr,2))/8.f + cpow(lr,2)*cexp(lr/2.) + 3.*lr*cexp(lr/2.)/2. - 1.)/cpow(lr,3);
	  b2 += h*(3.*lr - 2.*cexp(lr/2.0) + 1.25*cpow(lr,2) - 2.*lr*cexp(lr/2.) + 2.)/cpow(lr,3);
	  b3 +=-h*(lr - cexp(lr/2.) + 0.375*cpow(lr,2) - 0.5*lr*cexp(lr/2.) + 1.)/cpow(lr,3);
	}


	// Full dt coeeficients
	solver->MRSAAB_A[id+0] = creal(a1)/Nr;
	solver->MRSAAB_A[id+1] = creal(a2)/Nr;
	solver->MRSAAB_A[id+2] = creal(a3)/Nr;
	solver->MRSAAB_A[id+3] = 0.f;
	// Half coefficients
	solver->MRSAAB_B[id+0] = creal(b1)/Nr;
	solver->MRSAAB_B[id+1] = creal(b2)/Nr;
	solver->MRSAAB_B[id+2] = creal(b3)/Nr;
	solver->MRSAAB_B[id+3] = 0.f;

	// MRAB coefficients
	solver->MRAB_A[id+0]   =  23.*h/12. ;
	solver->MRAB_A[id+1]   = -16.*h/12. ;
	solver->MRAB_A[id+2]   =  5. *h/12. ;
	solver->MRAB_A[id+3]   =  0.f;

	solver->MRAB_B[id+0]   =  17.*h/24. ;
	solver->MRAB_B[id+1]   = - 7.*h/24. ;
	solver->MRAB_B[id+2]   =   2.*h/24. ;
	solver->MRAB_B[id+3]   =   0.f;
      }
      else {
	double complex a1 = 0. + 0.* I;
	double complex b1 = 0. + 0.* I;
	double complex a2 = 0. + 0.* I;
	double complex b2 = 0. + 0.* I;
	double complex a3 = 0. + 0.* I;
	double complex b3 = 0. + 0.* I;
	double complex a4 = 0. + 0.* I;
	double complex b4 = 0. + 0.* I;

	for(iint i = 0; i<Nr; ++i ){
	  double complex lr = alpha  + R[i];
	  a1 += h*(cexp(lr) - 3.*lr - (13.*cpow(lr,2))/3. - 4.*cpow(lr,3) + (11.*cpow(lr,2)*cexp(lr))/6. + cpow(lr,3)*cexp(lr) + 2.*lr*cexp(lr) - 1.)/cpow(lr,4);
	  a2 += h*(8.*lr - 3.*cexp(lr) + (19.*cpow(lr,2))/2. + 6.*cpow(lr,3) - 3.*cpow(lr,2)*cexp(lr) - 5.*lr*cexp(lr) + 3.)/cpow(lr,4);
	  a3 +=-h*(7.*lr - 3.*cexp(lr) + 7.*cpow(lr,2) + 4.*cpow(lr,3) - (3.*cpow(lr,2)*cexp(lr))/2. - 4.*lr*cexp(lr) + 3.)/cpow(lr,4);
	  a4 += h*(2.*lr - cexp(lr) + (11.*cpow(lr,2))/6. + cpow(lr,3) - (cpow(lr,2)*cexp(lr))/3. - lr*cexp(lr) + 1.)/cpow(lr,4);
	  b1 += h*(cexp((lr)/2.) - (5.*lr)/2. - (71.*cpow(lr,2))/24. - (35.*cpow(lr,3))/16. + (11.*cpow(lr,2)*cexp((lr)/2.))/6. + cpow(lr,3)*cexp((lr)/2.) + 2.*lr*cexp((lr)/2.) - 1.)/cpow(lr,4);
	  b2 += h*((13.*lr)/2. - 3.*cexp((lr)/2.) + (47.*cpow(lr,2))/8. + (35.*cpow(lr,3))/16. - 3.*cpow(lr,2)*cexp((lr)/2.) - 5.*lr*cexp((lr)/2.) + 3.)/cpow(lr,4);
	  b3 +=-h*((11.*lr)/2. - 3.*cexp((lr)/2.) + (31.*cpow(lr,2))/8. + (21.*cpow(lr,3))/16. - (3.*cpow(lr,2)*cexp((lr)/2.))/2. - 4.*lr*cexp((lr)/2.) + 3.)/cpow(lr,4);
	  b4 +=h*((3.*lr)/2. - cexp((lr)/2.) + (23.*cpow(lr,2))/24. + (5.*cpow(lr,3))/16. - (cpow(lr,2)*cexp((lr)/2.))/3. - lr*cexp((lr)/2.) + 1)/cpow(lr,4);
	    
	}


	// Full dt coeeficients
	solver->MRSAAB_A[id+0] = creal(a1)/Nr;
	solver->MRSAAB_A[id+1] = creal(a2)/Nr;
	solver->MRSAAB_A[id+2] = creal(a3)/Nr;
	solver->MRSAAB_A[id+3] = creal(a4)/Nr;
	// Half coefficients
	solver->MRSAAB_B[id+0] = creal(b1)/Nr;
	solver->MRSAAB_B[id+1] = creal(b2)/Nr;
	solver->MRSAAB_B[id+2] = creal(b3)/Nr;
	solver->MRSAAB_B[id+3] = creal(a4)/Nr;

	// MRAB coefficients
	solver->MRAB_A[id+0]   =  55.*h/24. ;
	solver->MRAB_A[id+1]   = -59.*h/24. ;
	solver->MRAB_A[id+2]   =  37. *h/24. ;
	solver->MRAB_A[id+3]   =  -3.*h/8;

	solver->MRAB_B[id+0]   =  99.*h/128. ;
	solver->MRAB_B[id+1]   = - 187.*h/384. ;
	solver->MRAB_B[id+2]   =   107.*h/384. ;
	solver->MRAB_B[id+3]   =   -25*h/384. ;
      }
    }

    // Exponential part
    solver->MRSAAB_C[l]    = exp(alpha);
    solver->MRAB_C[l]      =   h ;
  }
}

void advectionSetupMRSAABQuad3D (solver_t *solver) {
  
  mesh_t *mesh = solver->mesh;
  
    solver->Nrhs = 4;

    dfloat *q_zero = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*solver->Nfields,
				      sizeof(dfloat));
    solver->fQ = (dfloat*) calloc((mesh->totalHaloPairs+mesh->Nelements)*mesh->Np*solver->Nfields,
				  sizeof(dfloat));
    solver->rhsq = (dfloat*) calloc(mesh->Nelements*solver->Nrhs*mesh->Np*solver->Nfields,
				    sizeof(dfloat));

    rk4_coeffs(solver);

    dfloat cfl_small = 0.5; // depends on the stability region size (was .4, then 2)
    dfloat cfl_large = cfl_small;//((mesh->N)*cfl_small)/2;

    solver->resq = (dfloat*) calloc(mesh->Nelements*mesh->Np*solver->Nfields,
				    sizeof(dfloat));
    
    solver->localdt = (dfloat *) calloc(mesh->Nelements,sizeof(dfloat));
    
    dfloat glmin = 1e9, glmax = -1e9;
    // set time step
    for(iint e=0;e<mesh->Nelements;++e){
      dfloat lmin = 1e9, lmax = 0;
      for(iint f=0;f<mesh->Nfaces;++f){
	for(iint n=0;n<mesh->Nfp;++n){
	  iint sid = mesh->Nsgeo*mesh->Nfp*mesh->Nfaces*e + mesh->Nsgeo*mesh->Nfp*f+n;
	  
	  dfloat sJ   = mesh->sgeo[sid + mesh->Nq*SJID];
	  dfloat invJ = mesh->sgeo[sid + mesh->Nq*IJID];
	  
	  // A = 0.5*h*L
	  // => J*2 = 0.5*h*sJ*2
	  // => h = 2*J/sJ
	  
	  dfloat hest = 2./(sJ*invJ);
	  
	  lmin = mymin(lmin, hest);
	  lmax = mymax(lmax, hest);
	}
      }
      if (mesh->cubeDistance[e] == 0) {
	solver->localdt[e]  = cfl_small*lmin/((mesh->N+1.)*(mesh->N+1.));
      }
      else {
	solver->localdt[e]  = cfl_large*lmin/((mesh->N+1.)*(mesh->N+1.));
      }
      
      glmin = mymin(glmin, lmin);
      glmax = mymax(glmax, lmax);
      
    }

    //turn off for advection
    //dt = mymin(dt, cfl/mesh->tauInv);
    
    solver->finalTime = 5;
    solver->NtimeSteps = solver->finalTime/solver->dt;

    //remove this once mrab setup is integrated
    mesh->finalTime = solver->finalTime;
    
    iint maxLevels=100;
    meshMRABSetupQuad3D(mesh,solver->localdt,maxLevels);

    //remove this once mrab setup is integrated
    solver->dt = mesh->dt;
    solver->NtimeSteps = mesh->NtimeSteps;
    
    dfloat dt = solver->dt;
    
    solver->lev_updates = (iint *) calloc(mesh->MRABNlevels,sizeof(iint));
    solver->MRABshiftIndex = (iint*) calloc(mesh->MRABNlevels,sizeof(iint));
    solver->shift = (iint*) calloc(mesh->MRABNlevels,sizeof(iint));
    
    printf("cfl = %g %g\n", cfl_small,cfl_large);
    printf("dt = %g\n", dt);
    printf("max wave speed = %g\n", sqrt(3.)*solver->sqrtRT);
    printf("global h in [%g,%g]\n", glmin, glmax);
    
    // errorStep
    solver->errorStep = 10*mesh->Nq;
    
    printf("dt = %g\n", solver->dt);
    
    mrab4_coeffs(solver);

    occa::kernelInfo kernelInfo;
    
    advectionSetupOccaQuad3D(solver,&kernelInfo);
       
    solver->o_q =
      solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat), solver->q);

    solver->o_qpre =
      solver->device.malloc(mesh->Np*(mesh->totalHaloPairs+mesh->Nelements)*solver->Nfields*sizeof(dfloat), solver->q);

    solver->o_rhsq =
      solver->device.malloc(mesh->Np*solver->Nrhs*mesh->Nelements*solver->Nfields*sizeof(dfloat), solver->rhsq);
    
    solver->o_qPreCorr =
      solver->device.malloc(mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat),q_zero);
    
    solver->o_prerhsq =
      solver->device.malloc(mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat),q_zero);

    solver->o_qPreFilter =
      solver->device.malloc(mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat),q_zero);

    solver->o_qPreFiltered =
      solver->device.malloc(mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat),q_zero);

    solver->o_MRABlevels =
      solver->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*sizeof(iint),mesh->MRABlevel);
    
    solver->o_MRABelementIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    solver->o_MRABhaloIds = (occa::memory *) malloc(mesh->MRABNlevels*sizeof(occa::memory));
    solver->o_fQ =
      solver->device.malloc((mesh->Nelements+mesh->totalHaloPairs)*mesh->Np*solver->Nfields*sizeof(dfloat),solver->fQ);
    solver->o_lev_updates = solver->device.malloc(mesh->MRABNlevels*sizeof(iint),solver->lev_updates);
    solver->o_shift = solver->device.malloc(mesh->MRABNlevels*sizeof(iint),solver->MRABshiftIndex);
    solver->o_qFilter =
      solver->device.malloc(solver->Nrhs*mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat),solver->rhsq);
      
    solver->o_qFiltered =
      solver->device.malloc(solver->Nrhs*mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat),solver->rhsq);
    solver->o_qCorr =
      solver->device.malloc(solver->Nrhs*mesh->Nelements*solver->Nfields*mesh->Np*sizeof(dfloat),solver->rhsq);
    solver->o_resq =
      solver->device.malloc(mesh->Np*mesh->Nelements*solver->Nfields*sizeof(dfloat), solver->resq);

    for (iint lev = 0; lev < mesh->MRABNlevels;lev++) {
      if (mesh->MRABNelements[lev])
        solver->o_MRABelementIds[lev]
	  = solver->device.malloc(mesh->MRABNelements[lev]*sizeof(iint),
				mesh->MRABelementIds[lev]);
      if (mesh->MRABNhaloElements[lev])
	solver->o_MRABhaloIds[lev]
	  = solver->device.malloc(mesh->MRABNhaloElements[lev]*sizeof(iint),
				  mesh->MRABhaloIds[lev]);
    }

    solver->volumeKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/advectionVolumeQuad3D.okl",
					 "advectionVolumeMRSAABQuad3D",
					 kernelInfo);
    
    solver->volumePreKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/advectionVolumeQuad3D.okl",
					 "advectionVolumeLSERKQuad3D",
					 kernelInfo);
    solver->volumeCorrectionKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeCorrectionQuad3D.okl",
					 "boltzmannVolumeCorrectionMRSAABQuad3D",
					 kernelInfo);
    solver->volumeCorrPreKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannVolumeCorrectionQuad3D.okl",
					 "boltzmannVolumeCorrectionLSERKQuad3D",
					 kernelInfo);
    solver->surfaceKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/advectionSurfaceQuad3D.okl",
					 "advectionSurfaceMRSAABQuad3D",
					 kernelInfo);

    solver->surfacePreKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/advectionSurfaceQuad3D.okl",
					 "advectionSurfaceLSERKQuad3D",
					 kernelInfo);
  
    solver->traceUpdateKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
					   "boltzmannMRSAAB4TraceUpdateQuad3D",
					   kernelInfo);

    solver->traceUpdatePreKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
					 "boltzmannLSERKTraceUpdateQuad3D",
					 kernelInfo);

    solver->updateKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
					   "boltzmannMRSAAB4UpdateQuad3D",
					   kernelInfo);

    solver->updatePreKernel =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannUpdateQuad3D.okl",
					   "boltzmannLSERKUpdateQuad3D",
					   kernelInfo);
    
    solver->filterKernelq0H =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
					 "boltzmannFilterHq0Quad3D",
					 kernelInfo);
    solver->filterKernelq0V =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
					 "boltzmannFilterVq0Quad3D",
					 kernelInfo);

    solver->filterKernelH =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
					   "boltzmannFilterHQuad3D",
					   kernelInfo);
    solver->filterKernelV =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
					   "boltzmannFilterVQuad3D",
					   kernelInfo);
    
    solver->filterKernelHLSERK =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
					   "boltzmannFilterHLSERKQuad3D",
					   kernelInfo);
    solver->filterKernelVLSERK =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
					   "boltzmannFilterVLSERKQuad3D",
					   kernelInfo);
    
    solver->filterKernelLevelsH =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterHQuad3D.okl",
					   "boltzmannFilterLevelsHQuad3D",
					   kernelInfo);
    solver->filterKernelLevelsV =
      solver->device.buildKernelFromSource(DHOLMES "/okl/boltzmannFilterVQuad3D.okl",
					   "boltzmannFilterLevelsVQuad3D",
					   kernelInfo);
}
