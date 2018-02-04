#include "boltzmann2D.h"

// currently maximum
void boltzmannError2D(mesh2D *mesh, dfloat time, char *options){



 if(strstr(options,"PROBE")){
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      printf("Rank: %d has %d probes elements\n", rank, mesh->probeN);

      if(mesh->probeN){

      char fname[BUFSIZ];
      sprintf(fname, "ReynoldTest_%d_%.f_%05d_%04d_Alpha.dat", mesh->N, mesh->Re, mesh->Nelements, rank);

      FILE *fp; 
      fp = fopen(fname, "a");      
      
      fprintf(fp, "%2d %.4e %4e ", mesh->N, mesh->Re, time); 

      for(iint p=0; p<mesh->probeN; p++){
         iint e = mesh->probeElementIds[p];        
         dfloat srho = 0.0; 
         dfloat su = 0.0; 
         dfloat sv = 0.0; 
         for(iint n=0; n<mesh->Np; n++){
          dfloat rho  = mesh->q[mesh->Nfields*(n + e*mesh->Np) + 0];
          dfloat um   = mesh->q[mesh->Nfields*(n + e*mesh->Np) + 1]*mesh->sqrtRT/rho;
          dfloat vm   = mesh->q[mesh->Nfields*(n + e*mesh->Np) + 2]*mesh->sqrtRT/rho;
          // srho        += mesh->probeI[p*mesh->Np+n]*rho*mesh->RT;
          srho        += mesh->probeI[p*mesh->Np+n]*rho;
          su          += mesh->probeI[p*mesh->Np+n]*um;
          sv          += mesh->probeI[p*mesh->Np+n]*vm;
         }

        fprintf(fp, "%02d %.8e %.8e %.8e ",p, srho, su, sv); 
      }
    fprintf(fp, "\n"); 
    fclose(fp);
    }

  }







 if(strstr(options, "PML")){

    dfloat maxQ1 = 0, minQ1 = 1e9;
    iint fid = 0; //  
    for(iint e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        dfloat q1=0;
        iint id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];

        maxQ1 = mymax(maxQ1, fabs(mesh->q[id*mesh->Nfields + fid]));
        minQ1 = mymin(minQ1, fabs(mesh->q[id*mesh->Nfields + fid]));
      }
    }

    // compute maximum over all processes
    dfloat globalMaxQ1, globalMinQ1;
    MPI_Allreduce(&maxQ1, &globalMaxQ1, 1,
       MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minQ1, &globalMinQ1, 1,
      MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0)
      printf("%g %g %g (time,min(density),max(density)\n",
       time, globalMinQ1, globalMaxQ1);

    if(isnan(globalMinQ1) || isnan(globalMaxQ1))
      exit(EXIT_FAILURE);

}
else{

    // Coutte Flow exact solution for U velocity

    dfloat maxerr = 0, maxQ1 = 0, minQ1 = 1e9;
    iint fid = 1; 

    dfloat Uref        =  mesh->Ma*mesh->sqrtRT;
    dfloat nu = mesh->sqrtRT*mesh->sqrtRT/mesh->tauInv;

    for(iint e=0;e<mesh->Nelements;++e){
      for(iint n=0;n<mesh->Np;++n){
        dfloat q1=0;
        iint id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];
        // U = sqrt(RT)*Q2/Q1; 
       dfloat u   = mesh->sqrtRT*mesh->q[id*mesh->Nfields + 1]/mesh->q[id*mesh->Nfields+0];
 
        dfloat uex = y ; 
        for(iint k=1; k<=10; k++)
        {

         dfloat lamda = k*M_PI;
         // dfloat coef = -mesh->RT*mesh->tauInv/2. + sqrt(pow((mesh->RT*mesh->tauInv),2) /4.0 - (lamda*lamda*mesh->RT*mesh->RT));
         dfloat coef = -mesh->tauInv/2. + mesh->tauInv/2* sqrt(1.- 4.*pow(1./ mesh->tauInv, 2)* mesh->RT* lamda*lamda);
         uex += 2.*pow(-1,k)/(lamda)*exp(coef*time)*sin(lamda*y); //
         
         // uex += 2.*pow(-1,k)/(lamda)*exp(-nu*lamda*lamda*time)*sin(lamda*y); // !!!!!
        }

        maxerr = mymax(maxerr, fabs(u-uex));

        mesh->q[id*mesh->Nfields+2] = fabs(u-uex);

        maxQ1 = mymax(maxQ1, fabs(mesh->q[id*mesh->Nfields]));
        minQ1 = mymin(minQ1, fabs(mesh->q[id*mesh->Nfields]));
        
      }
    }

    // compute maximum over all processes
    dfloat globalMaxQ1, globalMinQ1, globalMaxErr;
    MPI_Allreduce(&maxQ1, &globalMaxQ1, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&minQ1, &globalMinQ1, 1, MPI_DFLOAT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&maxerr, &globalMaxErr, 1,  MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

    
    

   
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0){
      printf("%g %g %g %g (time,min(density),max(density),max(error)\n", time, globalMinQ1, globalMaxQ1, globalMaxErr);

      
      #if 1

    iint fld = 2;
    iint tstep = (time-mesh->startTime)/mesh->dt;
    char errname[BUFSIZ];
    sprintf(errname, "LSERK_err_long_%04d_%04d.vtu", rank, (tstep/mesh->errorStep));
    meshPlotVTU2D(mesh, errname,fld);
      
    iint tmethod = 0; 
    if(strstr(options,"LSERK"))
      tmethod = 1;
    if(strstr(options,"SRAB"))
       tmethod = 2;
    if(strstr(options,"SAAB"))
       tmethod = 3;
    if(strstr(options,"SARK"))
       tmethod = 4;
    if(strstr(options,"LSIMEX"))
       tmethod = 5;

      char fname[BUFSIZ]; 
      sprintf(fname, "boltzmannTemporalError.dat");
      FILE *fp; 
      fp = fopen(fname, "a");
      fprintf(fp, "%2d %2d %.8e %.8e %.8e %.8e %.8e\n", mesh->N, tmethod, mesh->dt, time,  globalMinQ1, globalMaxQ1, globalMaxErr); 
      fclose(fp); 

      mesh->maxErrorBoltzmann = globalMaxErr; 
      #endif
    }


    // if(isnan(globalMaxErr))
    //   exit(EXIT_FAILURE);
}

  
}
