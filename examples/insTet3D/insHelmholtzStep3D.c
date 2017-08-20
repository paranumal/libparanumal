#include "ins3D.h"

// complete a time step using LSERK4
void insHelmholtzStep3D(ins_t *ins, iint tstep,  iint haloBytes,
			dfloat * sendBuffer, dfloat * recvBuffer, 
			char   * options){
  
  mesh3D *mesh = ins->mesh; 
  
  solver_t *solver = ins->vSolver; 
  
  dfloat t = tstep*ins->dt + ins->dt;
  
  iint offset = mesh->Nelements+mesh->totalHaloPairs;

  iint rhsPackingMode = (strstr(options, "VECTORHELMHOLTZ")) ? 1:0;
  
  if(strstr(options,"SUBCYCLING")){
     // compute all forcing i.e. f^(n+1) - grad(Pr)
    ins->helmholtzRhsForcingKernel(mesh->Nelements,
				                          rhsPackingMode,
                                   mesh->o_vgeo,
                                   mesh->o_MM,
                                   ins->a0,
                                   ins->a1,
                                   ins->a2,
                                   ins->b0,
                                   ins->b1,
                                   ins->b2,
                                   ins->c0,
                                   ins->c1,
                                   ins->c2,
                                   ins->index,
                                   offset,
                                   ins->o_U,
                                   ins->o_V,
                                   ins->o_W,
                                   ins->o_NU,
                                   ins->o_NV,
                                   ins->o_NW,
                                   ins->o_Px,
                                   ins->o_Py,
                                   ins->o_Pz,
                                   ins->o_rhsU,
                                   ins->o_rhsV,
                                   ins->o_rhsW);
  }
  else{
       // compute all forcing i.e. f^(n+1) - grad(Pr)
    ins->helmholtzRhsForcingKernel(mesh->Nelements,
                                  rhsPackingMode,
                                   mesh->o_vgeo,
                                   mesh->o_MM,
                                   ins->a0,
                                   ins->a1,
                                   ins->a2,
                                   ins->b0,
                                   ins->b1,
                                   ins->b2,
                                   ins->c0,
                                   ins->c1,
                                   ins->c2,
                                   ins->index,
                                   offset,
                                   ins->o_U,
                                   ins->o_V,
                                   ins->o_W,
                                   ins->o_NU,
                                   ins->o_NV,
                                   ins->o_NW,
                                   ins->o_Px,
                                   ins->o_Py,
                                   ins->o_Pz,
                                   ins->o_rhsU,
                                   ins->o_rhsV,
                                   ins->o_rhsW);
  }
  
  ins->helmholtzRhsIpdgBCKernel(mesh->Nelements,
				                        rhsPackingMode,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                ins->tau,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                mesh->o_z,
                                mesh->o_vgeo,
                                mesh->o_sgeo,
                                mesh->o_EToB,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_DtT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                ins->o_rhsU,
                                ins->o_rhsV,
                                ins->o_rhsW);

  //use intermediate buffer for solve storage TODO: fix this later. Should be able to pull out proper buffer in elliptic solve
  if(rhsPackingMode==0){
    iint Ntotal = offset*mesh->Np;
    ins->o_UH.copyFrom(ins->o_U,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
    ins->o_VH.copyFrom(ins->o_V,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));
    ins->o_WH.copyFrom(ins->o_W,Ntotal*sizeof(dfloat),0,ins->index*Ntotal*sizeof(dfloat));

    #if 1
    printf("Solving for Ux \n");
    #if 1
      char fname[BUFSIZ];
      sprintf(fname, "/u0/outputs/ins3D/iterations.dat");
      FILE *fp;
      fp = fopen(fname, "a");
      fprintf(fp,"%d ",tstep);
      fclose(fp);
    #endif
    ellipticSolveTet3D( solver, ins->lambda, ins->o_rhsU, ins->o_UH, ins->vSolverOptions);
    
    printf("Solving for Uy \n");
    ellipticSolveTet3D(solver, ins->lambda, ins->o_rhsV, ins->o_VH, ins->vSolverOptions);
    
    printf("Solving for Uz \n");
    ellipticSolveTet3D(solver, ins->lambda, ins->o_rhsW, ins->o_WH, ins->vSolverOptions);
    #endif
    //copy into next stage's storage
    int index1 = (ins->index+1)%3; //hard coded for 3 stages
    ins->o_UH.copyTo(ins->o_U,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);
    ins->o_VH.copyTo(ins->o_V,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);  
    ins->o_WH.copyTo(ins->o_W,Ntotal*sizeof(dfloat),index1*Ntotal*sizeof(dfloat),0);  
  }else{

    printf("Solving for Ux,Uy and Uz \n");
    parAlmondPrecon(ins->o_rhsV, ins->precon->parAlmond, ins->o_rhsU); // rhs in rhsU, solution in rhsV

    iint Ntotal = mesh->Np*offset;
    dfloat *tmp  = (dfloat*) calloc(3*Ntotal, sizeof(dfloat));
    dfloat *tmpU = (dfloat*) calloc(Ntotal, sizeof(dfloat));
    dfloat *tmpV = (dfloat*) calloc(Ntotal, sizeof(dfloat));
    dfloat *tmpW = (dfloat*) calloc(Ntotal, sizeof(dfloat));
    
    ins->o_rhsV.copyTo(tmp);

    for(iint n=0;n<Ntotal;++n){
      tmpU[n] = tmp[3*n+0];
      tmpV[n] = tmp[3*n+1];
      tmpW[n] = tmp[3*n+2];
    }

    //copy into next stage's storage
    int index1 = (ins->index+1)%3; //hard coded for 3 stages
    
    ins->o_U.copyFrom(tmpU,Ntotal*sizeof(dfloat),index1*offset*mesh->Np*sizeof(dfloat));
    ins->o_V.copyFrom(tmpV,Ntotal*sizeof(dfloat),index1*offset*mesh->Np*sizeof(dfloat));
    ins->o_W.copyFrom(tmpW,Ntotal*sizeof(dfloat),index1*offset*mesh->Np*sizeof(dfloat));

    free(tmp);
    free(tmpU);
    free(tmpV);
    free(tmpW);
  }
  
}
