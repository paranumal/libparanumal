#include "ins3D.h"

// currently maximum
void insErrorNorms3D(ins_t *ins, dfloat time, char *options){

  mesh3D *mesh = ins->mesh;
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_V.copyTo(ins->V);  
  ins->o_W.copyTo(ins->W);  
  ins->o_P.copyTo(ins->P);

  dfloat maxU = 0.f;
  dfloat maxV = 0.f;
  dfloat maxW = 0.f;
  dfloat maxP = 0.f; 
  
  for(iint e=0;e<mesh->Nelements;++e){
    for(int n=0;n<mesh->Np;++n){
      iint id = n+e*mesh->Np;
      dfloat x = mesh->x[id];
      dfloat y = mesh->y[id];
      dfloat z = mesh->z[id];
      //
      id += ins->index*(mesh->Np)*(mesh->Nelements+mesh->totalHaloPairs);

      #if 0
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
      maxU = mymax(maxU, fabs(u-uex));
      maxV = mymax(maxV, fabs(v-vex));
      maxW = mymax(maxW, fabs(w-wex));
      maxP = mymax(maxP, fabs(p-pex));

      //
      ins->U[id] = fabs(u-uex);
      ins->V[id] = fabs(v-vex);
      ins->W[id] = fabs(w-wex);
      ins->P[id] = fabs(p-pex);
      // 
    }
  }





iint tstep = time/ins->dt+1; 
 int ranka;
 MPI_Comm_rank(MPI_COMM_WORLD, &ranka);
if(strstr(options, "VTU")){ 
    // output field files
    char fname[BUFSIZ];
    // sprintf(fname, "/u0/outputs/ins3D/");
    // sprintf(fname, "%sfoo_%04d", fname,rank);
    sprintf(fname, "/u0/outputs/ins3D/Error_%04d_%04d.vtu",ranka,tstep/ins->errorStep);
    
    insPlotVTU3D(ins, fname);
  } 



// compute maximum over all processes
dfloat gMaxU,gMaxV,gMaxW,gMaxP;
MPI_Allreduce(&maxU, &gMaxU, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
MPI_Allreduce(&maxV, &gMaxV, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
MPI_Allreduce(&maxW, &gMaxW, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);    
MPI_Allreduce(&maxP, &gMaxP, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);

int rank;
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
if(rank==0)
  printf("Step: %d Time: %g Uerr: %g Verr: %g Werr: %g Perr: %g\n", 
         (int)(time/ins->dt), time, gMaxU, gMaxV, gMaxW, gMaxP);

  // Do not Use mpi for Now!!!!!!!!!!!!!!!!!!!!!!1
  char fname[BUFSIZ];
  sprintf(fname, "/u0/outputs/ins3D/InfErr.dat");
  FILE *fp;
  fp = fopen(fname, "a");
  fprintf(fp,"Step: %d Time: %g Uerr: %g Verr: %g Werr: %g Perr: %g\n", 
         (int)(time/ins->dt), time, gMaxU, gMaxV, gMaxW, gMaxP);
 
  fclose(fp);
  




}
