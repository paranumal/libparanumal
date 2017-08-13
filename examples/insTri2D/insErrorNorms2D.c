#include "ins2D.h"

// currently maximum
void insErrorNorms2D(ins_t *ins, dfloat time, char *options){

  mesh2D *mesh = ins->mesh;
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V);  
  ins->o_P.copyTo(ins->P);


  #if 1

  const iint offset =  ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);
  //
  dfloat *dU = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *dV = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *dP = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  //
  dfloat *gU = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *gV = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));
  dfloat *gP = (dfloat*) calloc(mesh->Nelements*mesh->Np,sizeof(dfloat));

  dfloat nu = ins->nu;
  // Compute the Difference between exac solution and numerical one
  for(iint e=0;e<mesh->Nelements;++e){
    for(int i=0;i<mesh->Np;++i){
      iint id = i+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];

      // Compute Exact Solution
      #if 1
      dfloat uex = -sin(2.0 *M_PI*y)*exp(-nu*4.0*M_PI*M_PI*time);
      dfloat vex =  sin(2.0 *M_PI*x)*exp(-nu*4.0*M_PI*M_PI*time);
      dfloat pex = -cos(2.0 *M_PI*y)*cos(2.0*M_PI*x)*exp(-nu*8.0*M_PI*M_PI*time);
      //
      dfloat duexdx = 0.0;
      dfloat duexdy = -2.f*M_PI*cos(2.f*M_PI*y)*exp(-nu*4.0*M_PI*M_PI*time);; 
      dfloat dvexdx = 2.f*M_PI*cos(2.0*M_PI*x)*exp(-nu*4.0*M_PI*M_PI*time); ;
      dfloat dvexdy = 0.0;  
      dfloat dpexdx = 2.0*M_PI*cos(2.0 *M_PI*y)*sin(2.0*M_PI*x)*exp(-nu*8.0*M_PI*M_PI*time); 
      dfloat dpexdy = 2.0*M_PI*sin(2.0 *M_PI*y)*cos(2.0*M_PI*x)*exp(-nu*8.0*M_PI*M_PI*time); 
      #else
      dfloat lambda = 1./(2. * ins->nu) - sqrt(1./(4.*ins->nu * ins->nu) + 4.*M_PI*M_PI) ;
      dfloat uex = 1.0 - exp(lambda*x)*cos(2.*M_PI*y);
      dfloat vex =  lambda/(2.*M_PI)*exp(lambda*x)*sin(2.*M_PI*y);
      dfloat pex = 0.5*(1.0- exp(2.*lambda*x));
      //
      dfloat duexdx = -lambda*exp(lambda*x)*cos(2.0*M_PI*y);
      dfloat duexdy = 2.0*M_PI*exp(lambda*x)*sin(2.0*M_PI*y); 
      dfloat dvexdx = lambda*lambda/(2.0*M_PI)*exp(lambda*x)*sin(2.0*M_PI*y);
      dfloat dvexdy = lambda*exp(lambda*x)*cos(2.0*M_PI*y);  
      dfloat dpexdx = -lambda*exp(2.*lambda*x); 
      dfloat dpexdy = 0.0; 
      #endif

      // Compute the derivative of solution
      dfloat drdx = mesh->vgeo[e*mesh->Nvgeo + RXID];
      dfloat drdy = mesh->vgeo[e*mesh->Nvgeo + RYID];
      dfloat dsdx = mesh->vgeo[e*mesh->Nvgeo + SXID];
      dfloat dsdy = mesh->vgeo[e*mesh->Nvgeo + SYID];
      //
      dfloat dudr = 0, duds = 0;
      dfloat dvdr = 0, dvds = 0;
      dfloat dpdr = 0, dpds = 0; 

      for(int j=0;j<mesh->Np;++j){

        iint idj  = e*mesh->Np + j;

        dfloat u = ins->U[idj+offset];
        dfloat v = ins->V[idj+offset];
        dfloat p = ins->P[idj+offset];

        dfloat Drn = mesh->Dr[i*mesh->Np+j];
        dfloat Dsn = mesh->Ds[i*mesh->Np+j];

        dudr +=Drn*u; duds +=Dsn*u;
        dvdr +=Drn*v; dvds +=Dsn*v;
        dpdr +=Drn*p; dpds +=Dsn*p;
      }
      //  
      dfloat dudx = drdx*dudr + dsdx*duds;
      dfloat dudy = drdy*dudr + dsdy*duds;
      //
      dfloat dvdx = drdx*dvdr + dsdx*dvds;
      dfloat dvdy = drdy*dvdr + dsdy*dvds;
      //
      dfloat dpdx = drdx*dpdr + dsdx*dpds;
      dfloat dpdy = drdy*dpdr + dsdy*dpds;
      
      // Compute |grad(u_h)-grad(u_ex)|
      gU[id] = sqrt((dudx-duexdx)*(dudx-duexdx) + (dudy-duexdy)*(dudy-duexdy));
      gV[id] = sqrt((dvdx-dvexdx)*(dvdx-dvexdx) + (dvdy-dvexdy)*(dvdy-dvexdy));
      gP[id] = sqrt((dpdx-dpexdx)*(dpdx-dpexdx) + (dpdy-dpexdy)*(dpdy-dpexdy));
      //Store error field
      dU[id] = fabs(uex - ins->U[id+offset]);
      dV[id] = fabs(vex - ins->V[id+offset]);
      dP[id] = fabs(pex - ins->P[id+offset]);

    }
  }
  
  dfloat infu = 0,  infv= 0 , infp = 0; 
  dfloat l2u  = 0 , l2v = 0 , l2p  = 0; 
  dfloat h1u  = 0 , h1v = 0 , h1p  = 0; 
  for(iint e=0;e<mesh->Nelements;++e){
    dfloat l2ue=0, l2ve = 0, l2pe = 0;
    dfloat h1ue=0, h1ve = 0, h1pe = 0;

    for(int i=0;i<mesh->Np;++i){
      iint idi = i+e*mesh->Np;
      dfloat uei = dU[idi];
      dfloat vei = dV[idi];
      dfloat pei = dP[idi];
      //
      dfloat guei = gU[idi];
      dfloat gvei = gV[idi];
      dfloat gpei = gP[idi];

       for(int j=0;j<mesh->Np;++j){
        iint idj = j+e*mesh->Np;
        dfloat uej = dU[idj];
        dfloat vej = dV[idj];
        dfloat pej = dP[idj];
        //
        dfloat guej = gU[idj];
        dfloat gvej = gV[idj];
        dfloat gpej = gP[idj];

        dfloat mm = mesh->MM[j+i*mesh->Np];

        l2ue += mm*uei*uej;
        l2ve += mm*vei*vej;
        l2pe += mm*pei*pej;

        h1ue += mm*guei*guej;
        h1ve += mm*gvei*gvej;
        h1pe += mm*gpei*gpej;
       }
   
    // Infinity norm of error   
    infu = mymax(infu, uei);
    infv = mymax(infv, vei);
    infp = mymax(infp, pei);
    }

    dfloat j = mesh->vgeo[e*mesh->Nvgeo+JID];
    l2u += j*l2ue;
    l2v += j*l2ve;
    l2p += j*l2pe;
    //
    h1u += j*(h1ue + l2ue);
    h1v += j*(h1ve + l2ve);
    h1p += j*(h1pe + l2pe);
  }
  // Get Square Root
  l2u  = sqrt(l2u);
  l2v  = sqrt(l2v);
  l2p  = sqrt(l2p);
  //
  h1u  = sqrt(h1u);
  h1v  = sqrt(h1v);
  h1p  = sqrt(h1p);
  

   // compute maximum over all processes
    dfloat ginfu, ginfv, ginfp;
    dfloat gl2u, gl2v, gl2p;
    dfloat gh1u, gh1v, gh1p;

    MPI_Allreduce(&infu, &ginfu, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&infv, &ginfv, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&infp, &ginfp, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    //
    MPI_Allreduce(&l2u, &gl2u, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&l2v, &gl2v, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&l2p, &gl2p, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
     //
    MPI_Allreduce(&h1u, &gh1u, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&h1v, &gh1v, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&h1p, &gh1p, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);


     int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if(rank==0){
      char fname[BUFSIZ];
      sprintf(fname, "insErrors.txt");
      FILE *fp;
      fp = fopen(fname, "a");
      fprintf(fp,"%d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", 
                  mesh->N, ins->dt, ginfu, ginfv, ginfp,gl2u, gl2v, gl2p, gh1u, gh1v, gh1p);
      fclose(fp);
    }


 free(dU);
 free(dV);
 free(dP);
 free(gU);
 free(gV);
 free(gP);


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

// Do not Use mpi for Now!!!!!!!!!!!!!!!!!!!!!!1
char fname[BUFSIZ];
sprintf(fname, "/u0/outputs/ins2D/DragLift.dat");
FILE *fp;
fp = fopen(fname, "a");
fprintf(fp,"%d %.5e %.5e %.5e\n", mesh->N, time, cd, cl);
fclose(fp);

free(dUdx);
free(dVdx);
free(dUdy);
free(dVdy);

 #endif
  
}
