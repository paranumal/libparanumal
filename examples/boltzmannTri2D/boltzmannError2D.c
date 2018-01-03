#include "boltzmann2D.h"

// currently maximum
void boltzmannError2D(mesh2D *mesh, dfloat time,char *options){


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
    iint fid = 1; //U velocity

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
        for(iint k=1; k<=100; k++)
        {

         dfloat lamda = k*M_PI;
         uex += 2.*pow(-1,k)/(lamda)*exp(-nu*lamda*lamda*time)*sin(lamda*y);
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

    // iint fld = 2;
    // iint tstep = time/mesh->dt;
    // char errname[BUFSIZ];
    // sprintf(errname, "err_%04d_%04d.vtu", rank, (tstep/mesh->errorStep));
    // meshPlotVTU2D(mesh, errname,fld);
      
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
      fprintf(fp, "%2d %2d %.5e %.5e %.5e %.5e\n", mesh->N, tmethod, mesh->dt,  globalMinQ1, globalMaxQ1, globalMaxErr); 
      fclose(fp); 

      mesh->maxErrorBoltzmann = globalMaxErr; 
      #endif
    }


    // if(isnan(globalMaxErr))
    //   exit(EXIT_FAILURE);
}

  
}
