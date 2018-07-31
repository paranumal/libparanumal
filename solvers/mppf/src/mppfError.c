#include "mppf.h"


void mppfError(mppf_t *mppf, dfloat time){

  mesh_t *mesh = mppf->mesh;

  dlong offset = mesh->Np*(mesh->Nelements+mesh->totalHaloPairs);

  dfloat maxU   = 0, minU   = 1e9;
  dfloat maxV   = 0, minV   = 1e9;
  dfloat maxW   = 0, minW   = 1e9;
  dfloat maxP   = 0, minP   = 1e9; 
  dfloat maxR   = 0, minR   = 1e9; 
  dfloat maxM   = 0, minM   = 1e9; 
  dfloat maxPhi = 0, minPhi = 1e9; 

  if (mppf->options.compareArgs("EXACT", "NONE")) {// just compute maximums
 
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        int id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];
        dfloat z = mesh->z[id];

        maxU = mymax(maxU, fabs(mppf->U[id+0*offset]));
        minU = mymin(minU, fabs(mppf->U[id+0*offset]));
        
        maxV = mymax(maxV, fabs(mppf->U[id+1*offset]));
        minV = mymin(minV, fabs(mppf->U[id+1*offset]));
        
        if (mppf->dim==3) {
          maxW = mymax(maxW, fabs(mppf->U[id+3*offset]));
          minW = mymin(minW, fabs(mppf->U[id+3*offset]));  
        }

        maxP = mymax(maxP, fabs(mppf->P[id]));
        minP = mymin(minP, fabs(mppf->P[id]));

        maxR = mymax(maxR, fabs(mppf->Rho[id]));
        minR = mymin(minR, fabs(mppf->Rho[id]));

        maxM = mymax(maxM, fabs(mppf->Mu[id]));
        minM = mymin(minM, fabs(mppf->Mu[id]));

        maxPhi = mymax(maxPhi, fabs(mppf->Phi[id]));
        minPhi = mymin(minPhi, fabs(mppf->Phi[id]));

      }
    }

    // compute maximum over all processes
    dfloat gMaxU, gMinU;
    MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minU, &gMinU, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    dfloat gMaxV, gMinV;
    MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minV, &gMinV, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    dfloat gMaxW, gMinW;
    MPI_Allreduce(&maxW, &gMaxW, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minW, &gMinW, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
    
    dfloat gMaxP, gMinP;
    MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minP, &gMinP, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    dfloat gMaxR, gMinR;
    MPI_Allreduce(&maxR, &gMaxR, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minR, &gMinR, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    dfloat gMaxM, gMinM;
    MPI_Allreduce(&maxM, &gMaxM, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minM, &gMinM, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    dfloat gMaxPhi, gMinPhi;
    MPI_Allreduce(&maxPhi, &gMaxPhi, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minPhi, &gMinPhi, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);


    if(mesh->rank==0)
      if (mppf->dim==3) {
        printf("Step: %d Time: %g minU: %g maxU: %g minV: %g maxV: %g minW: %g maxW: %g minP: %g maxP: %g minRho: %g maxRho: %g minMu: %g maxMu: %g minPhi: %g maxPhi: %g \n", 
           (int)( (time-mppf->startTime)/mppf->dt)+1, time, gMinU, gMaxU, gMinV, gMaxV, gMinW, gMaxW, gMinP, gMaxP,
                  gMinR, gMaxR, gMinM, gMaxM, gMinPhi, gMaxPhi);
      } else {
        printf("Step: %d Time: %g minU: %g maxU: %g minV: %g maxV: %g minP: %g maxP: %g minRho: %g maxRho: %g minMu: %g maxMu: %g minPhi: %g maxPhi: %g \n", 
           (int)((time-mppf->startTime)/mppf->dt)+1, time, gMinU, gMaxU, gMinV, gMaxV, gMinP, gMaxP,
                  gMinR, gMaxR, gMinM, gMaxM, gMinPhi, gMaxPhi);
      }

    if( isnan(gMinU)   || isnan(gMaxU) || 
        isnan(gMinV)   || isnan(gMaxV) || 
        isnan(gMinW)   || isnan(gMaxW) || 
        isnan(gMinR)   || isnan(gMaxR) || 
        isnan(gMinM)   || isnan(gMaxM) || 
        isnan(gMinPhi) || isnan(gMaxPhi) ||
        isnan(gMinP)   || isnan(gMaxP) )
      exit(EXIT_FAILURE);
  } else { //compare to an exact solution

    if (mppf->options.compareArgs("EXACT","NSCH2D")) { //2D Taylor vortex

      for(dlong e=0;e<mesh->Nelements;++e){
        for(int n=0;n<mesh->Np;++n){
          int id = n+e*mesh->Np;
          dfloat x = mesh->x[id];
          dfloat y = mesh->y[id];

          dfloat uExact    =  cos(M_PI*y)*sin(M_PI*x)*sin(time);
          dfloat vExact    = -sin(M_PI*y)*cos(M_PI*x)*sin(time);
          dfloat pExact    =  sin(M_PI*y)*sin(M_PI*x)*cos(time);
          dfloat phiExact  =  cos(M_PI*x)*cos(M_PI*y)*sin(time);

          maxU   = mymax(maxU,   fabs(mppf->U[id+0*offset]-uExact));
          maxV   = mymax(maxV,   fabs(mppf->U[id+1*offset]-vExact));
          maxP   = mymax(maxP,   fabs(mppf->P[id]-pExact));
          maxPhi = mymax(maxPhi, fabs(mppf->Phi[id]-phiExact));

          #if 1
            mppf->U[id+0*offset] -= uExact;
            mppf->U[id+1*offset] -= vExact;
            mppf->P[id]          -= pExact;
            mppf->Phi[id]        -= phiExact;
          #endif
        }
      }
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      // compute maximum over all processes
      dfloat gMaxU;
      MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      dfloat gMaxV;
      MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
      
      dfloat gMaxP;
      MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      dfloat gMaxPhi;
      MPI_Allreduce(&maxPhi, &gMaxPhi, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

      if(mesh->rank==0)
        printf("Step: %d Time: %g ErrorU: %g ErrorV: %g ErrorP: %g ErrorPhi: %g \n", 
           (int)(time/mppf->dt), time, gMaxU, gMaxV, gMaxP, gMaxPhi);

      if( isnan(gMaxU) || 
          isnan(gMaxV) || 
          isnan(gMaxP) ||
          isnan(gMaxPhi))
        exit(EXIT_FAILURE);
    } 
    // else if (mppf->options.compareArgs("EXACT","BELTRAMI")) { //3D Beltrami flow

    //   for(dlong e=0;e<mesh->Nelements;++e){
    //     for(int n=0;n<mesh->Np;++n){
    //       dlong id = n+e*mesh->Np;
    //       dfloat x = mesh->x[id];
    //       dfloat y = mesh->y[id];
    //       dfloat z = mesh->z[id];

    //       dfloat a = M_PI/4.f;
    //       dfloat d = M_PI/2.f;

    //       dfloat uExact = -a*(exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y))*exp(-d*d*time);
    //       dfloat vExact = -a*(exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z))*exp(-d*d*time);
    //       dfloat wExact = -a*(exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x))*exp(-d*d*time);
    //       dfloat pExact = -a*a*exp(-2.f*d*d*time)*(exp(2.f*a*x)+exp(2.f*a*y)+exp(2.f*a*z))
    //                                               *( sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))
    //                                                 +sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(x+z))
    //                                                 +sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y))); 

    //       maxU = mymax(maxU, fabs(mppf->U[id+0*offset]-uExact));
    //       maxV = mymax(maxV, fabs(mppf->U[id+1*offset]-vExact));
    //       maxW = mymax(maxW, fabs(mppf->U[id+2*offset]-wExact));
    //       maxP = mymax(maxP, fabs(mppf->P[id]-pExact));

    //       #if 0
    //         mppf->U[id+0*offset] -= uExact;
    //         mppf->U[id+1*offset] -= vExact;
    //         mppf->U[id+2*offset] -= wExact;
    //         mppf->P[id] -= pExact;
    //       #endif
    //     }
    //   }

    //   // compute maximum over all processes
    //   dfloat gMaxU;
    //   MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

    //   dfloat gMaxV;
    //   MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

    //   dfloat gMaxW;
    //   MPI_Allreduce(&maxW, &gMaxW, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
      
    //   dfloat gMaxP;
    //   MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);

    //   if(mesh->rank==0)
    //     printf("Step: %d Time: %g ErrorU: %g ErrorV: %g ErrorW: %g ErrorP: %g \n", 
    //        (int)(time/mppf->dt), time, gMaxU, gMaxV, gMaxW, gMaxP);

    //   if( isnan(gMaxU) || 
    //       isnan(gMaxV) || 
    //       isnan(gMaxW) || 
    //       isnan(gMaxP) )
    //     exit(EXIT_FAILURE);
    // }
  }
}
