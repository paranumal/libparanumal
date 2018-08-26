/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "bns.h"

// TW: this needs to be tidied up

// currently maximum
void bnsError(bns_t *bns, dfloat time, setupAide &options){

 bns->o_q.copyTo(bns->q);

  mesh_t *mesh = bns->mesh; 

  // dfloat time = 0.0; 

  // if(options.compareArgs("TIME INTEGRATOR","MRSAAB"))
  //  time = bns->startTime + bns->dt*tstep*pow(2,(mesh->MRABNlevels-1));
  // else
  //  time = bns->startTime + tstep*bns->dt;

 if(bns->probeFlag){
    //   int rank;
    //   MPI_Comm_rank(mesh->comm, &rank);

    //   // printf("Rank: %d has %d probes elements\n", rank, mesh->probeN);

    //   if(mesh->probeN){

    //   char fname[BUFSIZ];
    //   sprintf(fname, "ProbeData_%d_%.f_%05d_%04d_%02d.dat", mesh->N, bns->Re, mesh->Nelements, rank, bns->Ntscale);

    //   FILE *fp; 
    //   fp = fopen(fname, "a");      
      
    //   fprintf(fp, "%2d %.4e %4e ", mesh->N, bns->Re, time); 

    //   for(int p=0; p<mesh->probeN; p++){
    //      int pid  = mesh->probeIds[p];
    //      int e    = mesh->probeElementIds[p];        
    //      dfloat srho = 0.0; 
    //      dfloat su = 0.0; 
    //      dfloat sv = 0.0; 
    //      for(int n=0; n<mesh->Np; n++){
    //       dfloat rho  = bns->q[bns->Nfields*(n + e*mesh->Np) + 0];
    //       dfloat um   = bns->q[bns->Nfields*(n + e*mesh->Np) + 1]*bns->sqrtRT/rho;
    //       dfloat vm   = bns->q[bns->Nfields*(n + e*mesh->Np) + 2]*bns->sqrtRT/rho;
    //       // srho        += mesh->probeI[p*mesh->Np+n]*rho*bns->RT;
    //       srho        += mesh->probeI[p*mesh->Np+n]*rho;
    //       su          += mesh->probeI[p*mesh->Np+n]*um;
    //       sv          += mesh->probeI[p*mesh->Np+n]*vm;
    //      }

    //     fprintf(fp, "%02d %04d %.8e %.8e %.8e ",pid, e, srho, su, sv); 
    //   }
    // fprintf(fp, "\n"); 
    // fclose(fp);
    // }

  }


  if(bns->outputForceStep)
    bnsForces(bns,time,options);



 if(options.compareArgs("ABSORBING LAYER","PML")){

    dfloat maxQ1 = 0, minQ1 = 1e9 , minQ3 = 1e9;
    int fid = 0; //  
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        dfloat q1=0, x=0. , y=0., z=0.;
        const dlong id = n+e*mesh->Np*bns->Nfields;
        maxQ1 = mymax(maxQ1, fabs(bns->q[id + fid*mesh->Np]));
        minQ1 = mymin(minQ1, fabs(bns->q[id + fid*mesh->Np]));
      }
    }

    // compute maximum over all processes
    dfloat globalMaxQ1, globalMinQ1;
    MPI_Allreduce(&maxQ1, &globalMaxQ1, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    MPI_Allreduce(&minQ1, &globalMinQ1, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);

    if(mesh->rank==0)
      printf("%g %g %g (time,min(density),max(density)\n", time, globalMinQ1, globalMaxQ1);

    if(isnan(globalMinQ1) || isnan(globalMaxQ1))
      exit(EXIT_FAILURE);

}


else{

    // // Coutte Flow exact solution for U velocity

    // dfloat maxerr = 0, maxQ1 = 0, minQ1 = 1e9;
    // int fid = 1; 

    // dfloat Uref        =  bns->Ma*bns->sqrtRT;
    // dfloat nu = bns->sqrtRT*bns->sqrtRT/bns->tauInv;

    // for(int e=0;e<mesh->Nelements;++e){
    //   for(int n=0;n<mesh->Np;++n){
    //     dfloat q1=0;
    //     int id = n+e*mesh->Np;
    //     dfloat x = mesh->x[id];
    //     dfloat y = mesh->y[id];
    //     // U = sqrt(RT)*Q2/Q1; 
    //    dfloat u   = bns->sqrtRT*bns->q[id*bns->Nfields + 1]/bns->q[id*bns->Nfields+0];
 
    //     dfloat uex = y ; 
    //     for(int k=1; k<=10; k++)
    //     {

    //      dfloat lamda = k*M_PI;
    //      // dfloat coef = -bns->RT*bns->tauInv/2. + sqrt(pow((bns->RT*bns->tauInv),2) /4.0 - (lamda*lamda*bns->RT*bns->RT));
    //      dfloat coef = -bns->tauInv/2. + bns->tauInv/2* sqrt(1.- 4.*pow(1./ bns->tauInv, 2)* bns->RT* lamda*lamda);
    //      uex += 2.*pow(-1,k)/(lamda)*exp(coef*time)*sin(lamda*y); //
         
    //      // uex += 2.*pow(-1,k)/(lamda)*exp(-nu*lamda*lamda*time)*sin(lamda*y); // !!!!!
    //     }

    //     maxerr = mymax(maxerr, fabs(u-uex));

    //     bns->q[id*bns->Nfields+2] = fabs(u-uex);

    //     maxQ1 = mymax(maxQ1, fabs(bns->q[id*bns->Nfields]));
    //     minQ1 = mymin(minQ1, fabs(bns->q[id*bns->Nfields]));
        
    //   }
    // }

    // // compute maximum over all processes
    // dfloat globalMaxQ1, globalMinQ1, globalMaxErr;
    // MPI_Allreduce(&maxQ1, &globalMaxQ1, 1, MPI_DFLOAT, MPI_MAX, mesh->comm);
    // MPI_Allreduce(&minQ1, &globalMinQ1, 1, MPI_DFLOAT, MPI_MIN, mesh->comm);
    // MPI_Allreduce(&maxerr, &globalMaxErr, 1,  MPI_DFLOAT, MPI_MAX, mesh->comm);

    
    

   
    // int rank;
    // MPI_Comm_rank(mesh->comm, &rank);
    // if(rank==0){
    //   printf("%g %g %g %g (time,min(density),max(density),max(error)\n", time, globalMinQ1, globalMaxQ1, globalMaxErr);
    // }


    // if(isnan(globalMaxErr))
    //   exit(EXIT_FAILURE);
}

  
}
