#include "insTet3D.h"

// currently maximum
void insErrorNormsTet3D(ins_t *ins, dfloat time, char *options){

  mesh3D *mesh = ins->mesh;
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V);  
  ins->o_W.copyTo(ins->W);  
  ins->o_P.copyTo(ins->P);

 
  dfloat *dU = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dV = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dW = (dfloat*) calloc(mesh->Np,sizeof(dfloat));
  dfloat *dP = (dfloat*) calloc(mesh->Np,sizeof(dfloat));

  dfloat l2u=0, l2v = 0, l2w =0, l2p = 0;
  dfloat liu=0, liv = 0, liw =0, lip = 0;

  for(dlong e=0;e<mesh->Nelements;++e){
    
    for(int n=0;n<mesh->Np;++n){
      dlong id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];
      //
      id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

      #if 1
      dfloat a = M_PI/4.0f, d = M_PI/2.0f; 
      dfloat uex = -a*( exp(a*x)*sin(a*y+d*z)+exp(a*z)*cos(a*x+d*y) )* exp(-d*d*time);
      dfloat vex = -a*( exp(a*y)*sin(a*z+d*x)+exp(a*x)*cos(a*y+d*z) )* exp(-d*d*time);
      dfloat wex = -a*( exp(a*z)*sin(a*x+d*y)+exp(a*y)*cos(a*z+d*x) )* exp(-d*d*time);

      dfloat pex= -a*a*exp(-2.f*d*d*time)*( exp(2.f*a*x) +exp(2.f*a*y)+exp(2.f*a*z))*( 
                        sin(a*x+d*y)*cos(a*z+d*x)*exp(a*(y+z))+
                        sin(a*y+d*z)*cos(a*x+d*y)*exp(a*(x+z))+
                        sin(a*z+d*x)*cos(a*y+d*z)*exp(a*(x+y))); 
      #else 
      dfloat lambda = 1./(2.*ins->nu)-sqrt(1./(4.*ins->nu*ins->nu) + 4.*M_PI*M_PI) ;
      dfloat a = M_PI/4.0f, d = M_PI/2.0f; 
      dfloat uex = 1.0 - exp(lambda*x)*cos(2.*M_PI*y);
      dfloat vex = lambda/(2.*M_PI)*exp(lambda*x)*sin(2.*M_PI*y);
      dfloat wex = 0.0;

      dfloat pex= 0.5*(1.0- exp(2.*lambda*x)); 
      #endif

      // Compute Maximum Nodal Error 
      dfloat u = ins->U[id];
      dfloat v = ins->V[id];
      dfloat w = ins->W[id];
      dfloat p = ins->P[id];
      // 
      liu = mymax(liu, fabs(u-uex));
      liv = mymax(liv, fabs(v-vex));
      liw = mymax(liw, fabs(w-wex));
      lip = mymax(lip, fabs(p-pex));
      //
      dU[n] = fabs(u-uex);
      dV[n] = fabs(v-vex);
      dW[n] = fabs(w-wex);
      dP[n] = fabs(p-pex);
      //
    }

    dfloat l2ue=0, l2ve = 0, l2we =0, l2pe = 0;
    
    for(int i=0;i<mesh->Np;++i){
    dfloat uei = dU[i];
    dfloat vei = dV[i];
    dfloat wei = dW[i];
    dfloat pei = dP[i];
      for(int j=0;j<mesh->Np;++j){
      dfloat uej = dU[j];
      dfloat vej = dV[j];
      dfloat wej = dW[j];
      dfloat pej = dP[j];
      dfloat mm = mesh->MM[j+i*mesh->Np];

      l2ue += mm*uei*uej;
      l2ve += mm*vei*vej;
      l2we += mm*wei*wej;
      l2pe += mm*pei*pej;
      }
    // 
    }

    dfloat j = mesh->vgeo[e*mesh->Nvgeo+JID];
    l2u += j*l2ue;
    l2v += j*l2ve;
    l2w += j*l2we;
    l2p += j*l2pe;
    //
  }
 
 
 free(dU);
 free(dV);
 free(dW);
 free(dP);


// compute maximum over all processes
dfloat gliu,gliv,gliw,glip;
MPI_Allreduce(&liu, &gliu, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
MPI_Allreduce(&liv, &gliv, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
MPI_Allreduce(&liw, &gliw, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);    
MPI_Allreduce(&lip, &glip, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
//
dfloat gl2u,gl2v,gl2w,gl2p;
MPI_Allreduce(&l2u, &gl2u, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&l2v, &gl2v, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&l2w, &gl2w, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
MPI_Allreduce(&l2p, &gl2p, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
//
// Get Square Root
  gl2u  = sqrt(gl2u);
  gl2v  = sqrt(gl2v);
  gl2w  = sqrt(gl2w);
  gl2p  = sqrt(gl2p);

#if 0
int rank;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if(rank==0){
  printf("Step: %d Time: %g Uerr: %g Verr: %g Werr: %g Perr: %g\n", 
         (int)(time/ins->dt), time, gliu, gliv, gliw, glip);

  char fname[BUFSIZ];
    // sprintf(fname, "insErrors.txt");
    //sprintf(fname, "beltrami_Ns%d.dat",ins->Nsubsteps);
    sprintf(fname, "BeltramiTemporal3D.txt");
    FILE *fp;
    fp = fopen(fname, "a");

    fprintf(fp,"%d %.5e %.5e %d %d %d %d %d %.5e %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", 
             mesh->N,ins->dt, time, ins->Nsubsteps, ins->NiterU, ins->NiterV, ins->NiterW, ins->NiterP, 
             gliu, gliv, gliw, glip, gl2u, gl2v, gl2w, gl2p);
    fclose(fp);  
 }
#endif 

}
