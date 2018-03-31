#include "advectionQuad3D.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=3){
    // to run cavity test case with degree N elements
    printf("usage: ./main meshes/cavityH005.msh N\n");
    exit(-1);
  }

  // int specify polynomial degree 
  int N = atoi(argv[2]);

  // set up mesh stuff
  dfloat sphereRadius = 1;
  mesh_t *mesh = meshSetupQuad3D(argv[1], N, sphereRadius);

  // set up boltzmann stuff
  solver_t *solver = advectionSetupMRQuad3D(mesh);

  /*  for(int e=0;e<mesh->Nelements;++e){
    for(int f=0;f<mesh->Nfaces;++f){
      for(int n=0;n<mesh->Nq;++n){
	int vidM = mesh->vmapM[e*mesh->Nfaces*mesh->Nq + f*mesh->Nq + n];
	int vidP = mesh->vmapP[e*mesh->Nfaces*mesh->Nq + f*mesh->Nq + n];
	if (mesh->x[vidM]-mesh->x[vidP] > 1e-12) printf("x mismatch = %lg\n", mesh->x[vidM]-mesh->x[vidP]);
	if (mesh->y[vidM]-mesh->y[vidP]) printf("y mismatch = %lg\n", mesh->y[vidM]-mesh->y[vidP]);
	if (mesh->z[vidM]-mesh->z[vidP]) printf("z mismatch = %lg\n", mesh->z[vidM]-mesh->z[vidP]);
      }
    }
    }*/
  
  // time step Boltzmann equations
  advectionRunLSERKQuad3D(solver);
  //advectionRunMRSAABQuad3D(solver);

  mesh->o_q.copyTo(mesh->q);
  dfloat l2 = 0;
  for (iint e = 0; e < mesh->Nelements; ++e) {
    for (iint n = 0; n < mesh->Np; ++n) {
      dfloat x = mesh->x[e*mesh->Np + n];
      dfloat y = mesh->y[e*mesh->Np + n];
      dfloat z = mesh->z[e*mesh->Np + n];
      dfloat t = mesh->finalTime;

      //rotate reference frame back to original
      dfloat xrot = x*cos(t) + y*sin(t);
      dfloat yrot = -1*x*sin(t) + y*cos(t);
      dfloat zrot = z;

      //current q0 is a gaussian pulse
      dfloat qref = 1 + .1*exp(-20*((xrot-1)*(xrot-1)+yrot*yrot+zrot*zrot));

      dfloat J = mesh->vgeo[mesh->Nvgeo*mesh->Np*e + JWID*mesh->Np + n];
      
      l2 += J*(qref - mesh->q[e*mesh->Np*mesh->Nfields + n])*(qref - mesh->q[e*mesh->Np*mesh->Nfields + n]);
    }
  }
  
  printf("norm %.5e\n",sqrt(l2));
  
  // close down MPI
  MPI_Finalize();

  exit(0);
  return 0;
}
