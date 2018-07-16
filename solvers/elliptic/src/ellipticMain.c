#include "elliptic.h"

int main(int argc, char **argv){

  // start up MPI
  MPI_Init(&argc, &argv);

  if(argc!=2){
    printf("usage: ./ellipticMain setupfile\n");

    MPI_Finalize();
    exit(-1);
  }

  // if argv > 2 then should load input data from argv
  setupAide options(argv[1]);
  
  // set up mesh stuff
  string fileName;
  int N, dim, elementType;

  options.getArgs("MESH FILE", fileName);
  options.getArgs("POLYNOMIAL DEGREE", N);
  options.getArgs("ELEMENT TYPE", elementType);
  options.getArgs("MESH DIMENSION", dim);

  // set up mesh
  mesh_t *mesh;
  switch(elementType){
  case TRIANGLES:
    mesh = meshSetupTri2D((char*)fileName.c_str(), N); break;
  case QUADRILATERALS:
    mesh = meshSetupQuad2D((char*)fileName.c_str(), N); break;
  case TETRAHEDRA:
    mesh = meshSetupTet3D((char*)fileName.c_str(), N); break;
  case HEXAHEDRA:
    mesh = meshSetupHex3D((char*)fileName.c_str(), N); break;
  }

  // parameter for elliptic problem (-laplacian + lambda)*q = f
  dfloat lambda;
  options.getArgs("LAMBDA", lambda);

  // set up
  occa::kernelInfo kernelInfo;
  elliptic_t *elliptic = ellipticSetup(mesh, lambda, kernelInfo, options);

  if(options.compareArgs("BENCHMARK", "BK5") ||
     options.compareArgs("BENCHMARK", "BP5")){
    
    // test Ax throughput
    occa::streamTag startAx = mesh->device.tagStream();
    
    int NAx = 1;
    
    for(int it=0;it<NAx;++it){
      // include gather-scatter
      if(options.compareArgs("BENCHMARK", "BP5"))
	ellipticOperator(elliptic, lambda, elliptic->o_x, elliptic->o_Ax, dfloatString); // standard precision

      if(options.compareArgs("BENCHMARK", "BK5")){
	if(!options.compareArgs("ELEMENT MAP", "TRILINEAR")){
	  elliptic->partialAxKernel(elliptic->NlocalGatherElements,			      
				    elliptic->o_localGatherElementList,
				    mesh->o_ggeo, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM,
				    lambda, elliptic->o_x, elliptic->o_Ax);
	}
	else{
	  elliptic->partialAxKernel(elliptic->NlocalGatherElements,			      
				    elliptic->o_localGatherElementList,
				    elliptic->o_EXYZ, elliptic->o_gllzw, mesh->o_Dmatrices, mesh->o_Smatrices, mesh->o_MM,
				    lambda, elliptic->o_x, elliptic->o_Ax);
	}
      }
    }
      
    occa::streamTag stopAx = mesh->device.tagStream();
      
    mesh->device.finish();
      
    double elapsedAx = mesh->device.timeBetween(startAx, stopAx);
    elapsedAx /= NAx;
      
      
    printf("%d, %d, %g, %d, %g, %g; \%\%elemental: N, dofs, elapsed, dummy, time per node, nodes/time %s\n",
	   mesh->N,
	   elliptic->NlocalGatherElements*mesh->Np,
	   0,
	   elapsedAx,
	   elapsedAx/(mesh->Np*mesh->Nelements),
	   mesh->Nelements*mesh->Np/elapsedAx,
	   options.getArgs("DISCRETIZATION").c_str());
      
  }
  else{
    
    // convergence tolerance
    dfloat tol = 1e-8;
  
    occa::streamTag startTag = mesh->device.tagStream();
  
    int it = ellipticSolve(elliptic, lambda, tol, elliptic->o_r, elliptic->o_x);

    occa::streamTag stopTag = mesh->device.tagStream();
    mesh->device.finish();
  
    double elapsed = mesh->device.timeBetween(startTag, stopTag);

    printf("%d, %d, %g, %d, %g, %g; \%\%global: N, dofs, elapsed, iterations, time per node, nodes*iterations/time %s\n",
	   mesh->N,
	   mesh->Nelements*mesh->Np,
	   elapsed,
	   it,
	   elapsed/(mesh->Np*mesh->Nelements),
	   mesh->Nelements*(it*mesh->Np/elapsed),
	   options.getArgs("PRECONDITIONER").c_str());

    if(options.compareArgs("DISCRETIZATION","CONTINUOUS")){
      dfloat zero = 0.;
      elliptic->addBCKernel(mesh->Nelements,
			    zero,
			    mesh->o_x,
			    mesh->o_y,
			    mesh->o_z,
			    elliptic->o_mapB,
			    elliptic->o_x);
    }
      
    // copy solution from DEVICE to HOST
    elliptic->o_x.copyTo(mesh->q);
      
    if (options.compareArgs("BASIS","BERN"))
      meshApplyElementMatrix(mesh,mesh->VB,mesh->q,mesh->q);
      
    dfloat maxError = 0;
    for(dlong e=0;e<mesh->Nelements;++e){
      for(int n=0;n<mesh->Np;++n){
	dlong   id = e*mesh->Np+n;
	dfloat xn = mesh->x[id];
	dfloat yn = mesh->y[id];
	dfloat zn = mesh->z[id];
      
	dfloat exact;
	if (elliptic->dim==2)
	  exact = sin(M_PI*xn)*sin(M_PI*yn);
	else 
	  exact = sin(M_PI*xn)*sin(M_PI*yn)*sin(M_PI*zn);
	dfloat error = fabs(exact-mesh->q[id]);
	  
	maxError = mymax(maxError, error);
	//mesh->q[id] -= exact;
      }
    }
      
    dfloat globalMaxError = 0;
    MPI_Allreduce(&maxError, &globalMaxError, 1, MPI_DFLOAT, MPI_MAX, MPI_COMM_WORLD);
    if(mesh->rank==0)
      fprintf(stderr,"globalMaxError = %g\n", globalMaxError);
      
#if 0
    char fname[BUFSIZ];
    string outName;
    options.getArgs("OUTPUT FILE NAME", outName);
    sprintf(fname, "%s_%04d.vtu",(char*)outName.c_str(), rank);
    if(elliptic->dim==3)
      meshPlotVTU3D(mesh, fname, 0);
    else 
      meshPlotVTU2D(mesh, fname, 0);
#endif
  }
  
  
  // close down MPI
  MPI_Finalize();
  
  exit(0);
  return 0;
}
