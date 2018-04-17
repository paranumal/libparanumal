#include "boltzmann2D.h"


void boltzmannForces2D(bns_t *bns, dfloat time, char * options){


  mesh2D *mesh = bns->mesh; 
	
	// Hard coded weights i.e. sum(MM1D,1), MM1D = inv(V1)' * inv(V1)
  dfloat W[mesh->Nfp];
  if(mesh->N==1){
    W[0] = 1.000000000000000; W[1] = 1.000000000000000; 
  }

 if(mesh->N==2){
    W[0] = 0.333333333333333; W[1] = 1.333333333333333;
    W[2] = 0.333333333333333; 
  }

  if(mesh->N==3){
    W[0] = 0.166666666666666; W[1] = 0.833333333333334;
    W[2] = 0.833333333333334; W[3] = 0.166666666666666;
  }
  if(mesh->N==4){
    W[0] = 0.100000000000000; W[1] = 0.544444444444444;
    W[2] = 0.711111111111111; W[3] = 0.544444444444444;
    W[4] = 0.100000000000000; 
  }
  if(mesh->N==5){
    W[0] = 0.066666666666666; W[1] = 0.378474956297847;
    W[2] = 0.554858377035487; W[3] = 0.554858377035486;
    W[4] = 0.378474956297847; W[5] = 0.066666666666667;
  }

  if(mesh->N==6){
    W[0] = 0.047619047619047; W[1] = 0.276826047361566;
    W[2] = 0.431745381209863; W[3] = 0.487619047619048;
    W[4] = 0.431745381209863; W[5] = 0.276826047361566;
    W[6] = 0.047619047619047;
  }

   

   dfloat Fx =0. , Fy = 0. ;
  //
   for(int e=0;e<mesh->Nelements;++e){
     int flag = 0; 
    for(int f=0;f<mesh->Nfaces;++f){
      int bc = mesh->EToB[e*mesh->Nfaces+f];
      if(bc == 1){ flag = 1; }
    }

    if(flag){

			for(int f=0;f<mesh->Nfaces;++f){
				int bc = mesh->EToB[e*mesh->Nfaces+f];
				if(bc==1){
					for(int n=0;n<mesh->Nfp; n++){
						int sid  = mesh->Nsgeo*(e*mesh->Nfaces+f);
						dfloat nx = mesh->sgeo[sid+NXID];
						dfloat ny = mesh->sgeo[sid+NYID];
						dfloat sJ = mesh->sgeo[sid+SJID];
						//
						int vid  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
            int idM  = mesh->vmapM[vid];
		    	 
						dfloat q1  = bns->q[mesh->Nfields*idM + 0];
						dfloat q2  = bns->q[mesh->Nfields*idM + 1];
						dfloat q3  = bns->q[mesh->Nfields*idM + 2];
						dfloat q4  = bns->q[mesh->Nfields*idM + 3];
						dfloat q5  = bns->q[mesh->Nfields*idM + 4];
						dfloat q6  = bns->q[mesh->Nfields*idM + 5];

						
    	      // Compute Stress Tensor
    	      dfloat s11 = -bns->RT*(sqrt(2.0)*q5 - q2*q2/q1);
    	      dfloat s12 = -bns->RT*(          q4 - q2*q3/q1);
    	      dfloat s22 = -bns->RT*(sqrt(2.0)*q6 - q3*q3/q1);

    	      dfloat P   = q1*bns->RT; // rho*RT 

            Fx  += W[n]*sJ*(P*nx  - (s11*nx + s12*ny) );
            Fy  += W[n]*sJ*(P*ny  - (s12*nx + s22*ny) );
          }
        }
      }
    }
  }


  dfloat gFx = 0., gFy = 0.; 
  // Add all processors force, 
  MPI_Allreduce(&Fx, &gFx, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Fy, &gFy, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  //
  int rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank==0){
		char fname[BUFSIZ];
		sprintf(fname, "ForceData_%d_%.f_%05d_%03d.dat", mesh->N, bns->Re, mesh->Nelements, bns->Ntscale);

		FILE *fp; 
		fp = fopen(fname, "a");  

		fprintf(fp, "%.4e %.8e %.8e \n", time, gFx, gFy); 

		fclose(fp);
  }
  
  




}