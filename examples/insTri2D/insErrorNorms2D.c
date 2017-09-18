#include "ins2D.h"

// currently maximum
void insErrorNorms2D(ins_t *ins, dfloat time, char *options){

  mesh2D *mesh = ins->mesh;
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V);  
  ins->o_P.copyTo(ins->P);


  #if 1

    dfloat *dU = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
    dfloat *dV = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
    dfloat *dP = (dfloat*) calloc(mesh->Np,sizeof(dfloat));

    dfloat l2u=0, l2v = 0, l2p = 0;
    dfloat liu=0, liv = 0, lip = 0;

    dfloat nu = ins->nu;
  
  for(iint e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      iint id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

      // Compute Exact Solution
      #if 1
      dfloat uex = -sin(2.0 *M_PI*y)*exp(-nu*4.0*M_PI*M_PI*time);
      dfloat vex =  sin(2.0 *M_PI*x)*exp(-nu*4.0*M_PI*M_PI*time);
      dfloat pex = -cos(2.0 *M_PI*y)*cos(2.0*M_PI*x)*exp(-nu*8.0*M_PI*M_PI*time);
     
      #else
      dfloat lambda = 1./(2. * ins->nu) - sqrt(1./(4.*ins->nu * ins->nu) + 4.*M_PI*M_PI) ;
      dfloat uex = 1.0 - exp(lambda*x)*cos(2.*M_PI*y);
      dfloat vex =  lambda/(2.*M_PI)*exp(lambda*x)*sin(2.*M_PI*y);
      dfloat pex = 0.5*(1.0- exp(2.*lambda*x));
      #endif
     // Compute Maximum Nodal Error 
      dfloat u = ins->U[id];
      dfloat v = ins->V[id];
      dfloat p = ins->P[id];
      // 
      liu = mymax(liu, fabs(u-uex));
      liv = mymax(liv, fabs(v-vex));
      lip = mymax(lip, fabs(p-pex));
      //
      dU[n] = fabs(u-uex);
      dV[n] = fabs(v-vex);
      dP[n] = fabs(p-pex);
      //
    }

    dfloat l2ue=0, l2ve = 0,  l2pe = 0;

    for(int i=0;i<mesh->Np;++i){
      dfloat uei = dU[i];
      dfloat vei = dV[i];
      dfloat pei = dP[i];

      for(int j=0;j<mesh->Np;++j){
        dfloat uej = dU[j];
        dfloat vej = dV[j];
        dfloat pej = dP[j];
        dfloat mm = mesh->MM[j+i*mesh->Np];

        l2ue += mm*uei*uej;
        l2ve += mm*vei*vej;
        l2pe += mm*pei*pej;
      }
    // 
    }

    dfloat j = mesh->vgeo[e*mesh->Nvgeo+JID];
    l2u += j*l2ue;
    l2v += j*l2ve;
    l2p += j*l2pe;
    //
  }

  free(dU);
  free(dV);
  free(dP);
      

  // compute maximum over all processes
  dfloat giu  = 0, giv = 0 , gip = 0;
  dfloat glu  = 0, glv = 0,  glp  = 0;

  MPI_Allreduce(&liu, &giu, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&liv, &giv, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&lip, &gip, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
  //
  MPI_Allreduce(&l2u, &glu, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&l2v, &glv, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&l2p, &glp, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);

  glu = sqrt(glu);
  glv = sqrt(glv);
  glp = sqrt(glp);


  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if(rank==0){
    printf("%d %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", mesh->N, ins->dt, giu, giv, gip, glu, glv, glp);

    char fname[BUFSIZ];
    // sprintf(fname, "insErrors.txt");
    sprintf(fname, "vortex_pamg_Ns%d_N%d.dat",ins->Nsubsteps, mesh->N);
    FILE *fp;
    fp = fopen(fname, "a");

    fprintf(fp,"%d %.5e %d %d %d %d %.5e %.5e %.5e %.5e %.5e %.5e\n", 
             mesh->N, time, ins->Nsubsteps, ins->NiterU, ins->NiterV, ins->NiterP, giu, giv, gip, glu, glv, glp);
    fclose(fp);

  }

 #else

 const iint offset =  ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);
  //
  dfloat *dUdx = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *dUdy = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *dVdx = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *dVdy = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  //
   for(iint e=0;e<mesh->Nelements;++e){
    //  int flag = 0; 
    // for(iint f=0;f<mesh->Nfaces;++f){
    //   iint bc = mesh->EToB[e*mesh->Nfaces+f];
    //   if(bc == 1){ flag = 1; }
    // }

    // if(flag){ // Compute Derivatives dUdx, etc;
      // prefetch geometric factors (constant on triangle)
      dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
      dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
      dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
      dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];
      for(iint n=0;n<mesh->Np;++n){
        dfloat dudr = 0, duds = 0, dvdr = 0, dvds = 0;     
        for(iint i=0;i<mesh->Np;++i){
          // load data at node i of element e (note Nfields==4)
          iint id = e*mesh->Np + i;
          dfloat u = ins->U[id+offset];
          dfloat v = ins->V[id+offset];
             
          dfloat Drni = mesh->Dr[n*mesh->Np+i];
          dfloat Dsni = mesh->Ds[n*mesh->Np+i];
          // differentiate (u,v,p) with respect to 'r' and 's'
          dudr += Drni*u;
          duds += Dsni*u;
          dvdr += Drni*v;
          dvds += Dsni*v;   
      }
      dUdx[e*mesh->Np + n] = drdx*dudr + dsdx*duds;
      dUdy[e*mesh->Np + n] = drdy*dudr + dsdy*duds;
      dVdx[e*mesh->Np + n] = drdx*dvdr + dsdx*dvds;
      dVdy[e*mesh->Np + n] = drdy*dvdr + dsdy*dvds;
    // }
  }

 }
 

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
//
dfloat cd = 0.0, cl = 0.0; 
iint sk=0;
for(iint e=0;e<mesh->Nelements;++e){
   iint flag = 0; 
   iint bc = 0; 
  for(iint f=0;f<mesh->Nfaces;++f){
    bc = mesh->EToB[e*mesh->Nfaces+f];
    if(bc == 1){ flag = 1; }
  }

  if(flag){
    for(iint f=0;f<mesh->Nfaces;++f){
      bc = mesh->EToB[e*mesh->Nfaces+f];
      if(bc==1){
       for(iint n=0;n<mesh->Nfp; n++){
        // load surface geofactors for this face
        iint sid = mesh->Nsgeo*(e*mesh->Nfaces+f);
        dfloat nx = mesh->sgeo[sid+0];
        dfloat ny = mesh->sgeo[sid+1];
        dfloat sJ = mesh->sgeo[sid+2];
        
        iint vid  = e*mesh->Nfp*mesh->Nfaces + f*mesh->Nfp + n;
        iint idM = mesh->vmapM[vid];
        //
        dfloat dudx = dUdx[idM]; 
        dfloat dudy = dUdy[idM]; 
        dfloat dvdx = dVdx[idM]; 
        dfloat dvdy = dVdy[idM];
        //
        dfloat p = ins->P[idM + offset];

        cd += W[n]*(sJ*(-p*nx + ins->nu*(nx*2.0*dudx + ny*(dvdx + dudy))));
        cl += W[n]*(sJ*(-p*ny + ins->nu*(nx*(dvdx + dudy) + ny*2.0*dvdy)));
      }
      // printf("%d %.5e \n",sk,cl);
      //sk++;
    }
  }
}
}

int rank;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if(rank==0){
  // Do not Use mpi for Now!!!!!!!!!!!!!!!!!!!!!!1
  char fname[BUFSIZ];
  //sprintf(fname, "/u0/outputs/ins2D/Report.dat");
  sprintf(fname, "report_Ns%d_N%d.dat",ins->Nsubsteps, mesh->N);


  FILE *fp;
  fp = fopen(fname, "a");
  fprintf(fp,"%d %.5e %d %d %d %d %.5e %.5e\n", mesh->N, time, ins->Nsubsteps, ins->NiterU, ins->NiterV, ins->NiterP, cd, cl);
  fclose(fp);
}

free(dUdx);
free(dVdx);
free(dUdy);
free(dVdy);

 #endif
  
}
