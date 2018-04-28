#include "boltzmann2D.h"

// currently maximum
void boltzmannError2D(bns_t *bns, dfloat time, char *options){

  mesh2D *mesh = bns->mesh; 


#if 1
 if(strstr(options,"PROBE")){
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      // printf("Rank: %d has %d probes elements\n", rank, mesh->probeN);

      if(mesh->probeN){

      char fname[BUFSIZ];
      sprintf(fname, "ProbeData_%d_%.f_%05d_%04d_%02d.dat", mesh->N, bns->Re, mesh->Nelements, rank, bns->Ntscale);

      FILE *fp; 
      fp = fopen(fname, "a");      
      
      fprintf(fp, "%2d %.4e %4e ", mesh->N, bns->Re, time); 

      for(int p=0; p<mesh->probeN; p++){
         int pid  = mesh->probeIds[p];
         int e    = mesh->probeElementIds[p];        
         dfloat srho = 0.0; 
         dfloat su = 0.0; 
         dfloat sv = 0.0; 
         for(int n=0; n<mesh->Np; n++){
          dfloat rho  = bns->q[bns->Nfields*(n + e*mesh->Np) + 0];
          dfloat um   = bns->q[bns->Nfields*(n + e*mesh->Np) + 1]*bns->sqrtRT/rho;
          dfloat vm   = bns->q[bns->Nfields*(n + e*mesh->Np) + 2]*bns->sqrtRT/rho;
          // srho        += mesh->probeI[p*mesh->Np+n]*rho*bns->RT;
          srho        += mesh->probeI[p*mesh->Np+n]*rho;
          su          += mesh->probeI[p*mesh->Np+n]*um;
          sv          += mesh->probeI[p*mesh->Np+n]*vm;
         }

        fprintf(fp, "%02d %04d %.8e %.8e %.8e ",pid, e, srho, su, sv); 
      }
    fprintf(fp, "\n"); 
    fclose(fp);
    }

  }








#else 
  if(strstr(options,"PROBE")){
    
    // Move this routine to probe setup



    int rank, size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int root = 0; // root rank
    
    int probeNfields = 3; // rho, u, v
    
    dfloat *probeData   = (dfloat *) malloc(mesh->probeN*probeNfields*sizeof(dfloat)); 
    dfloat *recvData    = (dfloat *) malloc(mesh->probeNTotal*probeNfields*sizeof(dfloat));
    int   *recvcount   = (int   *) malloc(size*sizeof(int));
    int   *recvdisp    = (int   *) malloc(size*sizeof(int));
    int   *probeIds    = (int   *) malloc(mesh->probeNTotal*sizeof(int)); 
    
    // Collect number of probes on the root  
    MPI_Gather(&(mesh->probeN), 1, MPI_int, recvcount, 1, MPI_int, root, MPI_COMM_WORLD);
    
    if(rank==root){
    recvdisp[0]= 0;  
      for(int i=1; i<size; ++i){
            recvdisp[i]   = recvdisp[i-1] + recvcount[i-1]; 
          //printf("%d %d \n", recvcount[i], recvdisp[i]);
      }
    }

    // Collect probe ids
    MPI_Gatherv(mesh->probeIds, mesh->probeN, MPI_int, probeIds, recvcount, recvdisp, MPI_int, root, MPI_COMM_WORLD);

    // if(rank==0){
    //   for(int id=0; id<mesh->probeNTotal; id++){
    //     printf("id: %d  ProbeId: %d \n", id, probeIds[id]);
    //   }
    // }

    // // Sort Probe Ids such that probes are ardered in ascending order

    // int *probeIndex = (int *) malloc(2*mesh->probeNTotal*sizeof(int));

    // for(int id=0; id<mesh->probeNTotal; id++){
    // probeIndex



    // }

    
    if(rank==root){
    recvdisp[0]= 0;  
      for(int i=0; i<size; ++i){
        recvcount[i] *=probeNfields;
          if(i>0)
            recvdisp[i]   = recvdisp[i-1] + recvcount[i-1]; 
          printf("%d %d \n", recvcount[i], recvdisp[i]);
      }
    }
    
  
    // fill probe Data
    if(mesh->probeN){
      for(int p=0; p<mesh->probeN; p++){
        int pid  = mesh->probeIds[p];
        int e    = mesh->probeElementIds[p];        
        dfloat sr = 0.0; 
        dfloat su = 0.0; 
        dfloat sv = 0.0; 
        for(int n=0; n<mesh->Np; n++){
          dfloat rho  = bns->q[bns->Nfields*(n + e*mesh->Np) + 0];
          dfloat um   = bns->q[bns->Nfields*(n + e*mesh->Np) + 1]*bns->sqrtRT/rho;
          dfloat vm   = bns->q[bns->Nfields*(n + e*mesh->Np) + 2]*bns->sqrtRT/rho;
          //
          sr += mesh->probeI[p*mesh->Np+n]*rho;
          su += mesh->probeI[p*mesh->Np+n]*um;
          sv += mesh->probeI[p*mesh->Np+n]*vm;
        }
        probeData[p*probeNfields + 0] = sr;
        probeData[p*probeNfields + 1] = su;
        probeData[p*probeNfields + 2] = sv;
      }
    }



    
    MPI_Gatherv(probeData, mesh->probeN*probeNfields, MPI_DFLOAT, recvData, recvcount, recvdisp, MPI_DFLOAT, root, MPI_COMM_WORLD);

    if(rank==root){     
       char fname[BUFSIZ];
      sprintf(fname, "ProbeData_%d_%.f_%05d.dat", mesh->N, mesh->Re, mesh->Nelements);

      FILE *fp; 
      fp = fopen(fname, "a");      

      fprintf(fp, "%4e ", time); 
      for(int p=0; p<mesh->probeNTotal; p++){
        fprintf(fp, "%02d %.8e %.8e %.8e ", probeIds[p],recvData[p*probeNfields+0], 
                                                        recvData[p*probeNfields+1],
                                                        recvData[p*probeNfields+2]);
      }
      fprintf(fp, "\n"); 
      fclose(fp);
    }
  }



#endif




 if(strstr(options, "PML")){

    dfloat maxQ1 = 0, minQ1 = 1e9;
    int fid = 0; //  
    for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        dfloat q1=0;
        int id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];

        maxQ1 = mymax(maxQ1, fabs(bns->q[id*bns->Nfields + fid]));
        minQ1 = mymin(minQ1, fabs(bns->q[id*bns->Nfields + fid]));
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
    int fid = 1; 

    dfloat Uref        =  bns->Ma*bns->sqrtRT;
    dfloat nu = bns->sqrtRT*bns->sqrtRT/bns->tauInv;

    for(int e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
        dfloat q1=0;
        int id = n+e*mesh->Np;
        dfloat x = mesh->x[id];
        dfloat y = mesh->y[id];
        // U = sqrt(RT)*Q2/Q1; 
       dfloat u   = bns->sqrtRT*bns->q[id*bns->Nfields + 1]/bns->q[id*bns->Nfields+0];
 
        dfloat uex = y ; 
        for(int k=1; k<=10; k++)
        {

         dfloat lamda = k*M_PI;
         // dfloat coef = -bns->RT*bns->tauInv/2. + sqrt(pow((bns->RT*bns->tauInv),2) /4.0 - (lamda*lamda*bns->RT*bns->RT));
         dfloat coef = -bns->tauInv/2. + bns->tauInv/2* sqrt(1.- 4.*pow(1./ bns->tauInv, 2)* bns->RT* lamda*lamda);
         uex += 2.*pow(-1,k)/(lamda)*exp(coef*time)*sin(lamda*y); //
         
         // uex += 2.*pow(-1,k)/(lamda)*exp(-nu*lamda*lamda*time)*sin(lamda*y); // !!!!!
        }

        maxerr = mymax(maxerr, fabs(u-uex));

        bns->q[id*bns->Nfields+2] = fabs(u-uex);

        maxQ1 = mymax(maxQ1, fabs(bns->q[id*bns->Nfields]));
        minQ1 = mymin(minQ1, fabs(bns->q[id*bns->Nfields]));
        
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

    int fld = 2;
    int tstep = (time-bns->startTime)/bns->dt;
    char errname[BUFSIZ];
    sprintf(errname, "LSERK_err_long_%04d_%04d.vtu", rank, (tstep/mesh->errorStep));
    meshPlotVTU2D(mesh, errname,fld);
      
    int tmethod = 0; 
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
      fprintf(fp, "%2d %2d %.8e %.8e %.8e %.8e %.8e\n", mesh->N, tmethod, bns->dt, time,  globalMinQ1, globalMaxQ1, globalMaxErr); 
      fclose(fp); 

      bns->maxError = globalMaxErr; 
      #endif
    }


    // if(isnan(globalMaxErr))
    //   exit(EXIT_FAILURE);
}

  
}
