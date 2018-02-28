#include "ins2D.h"

// complete a time step using LSERK4
void insPoissonStep2D(ins_t *ins, int tstep, int haloBytes,
				       dfloat * sendBuffer, dfloat * recvBuffer,
				        char   * options){

  mesh2D *mesh = ins->mesh;
  solver_t *solver = ins->pSolver;
  dfloat t = tstep*ins->dt + ins->dt;

  //hard coded for 3 stages.
  //The result of the helmholtz solve is stored in the next index
  int index1   = (ins->index+1)%3;
  int offset  = mesh->Nelements+mesh->totalHaloPairs;
  int ioffset = index1*offset;

  /* note: the surface kernel isn't needed with continuous pressure. Just the inflow boundary 
           contributions to the surface 
           TODO: Need a separate kernel to do the surface kernel for just boundaries */
  //if (strstr(ins->pSolverOptions,"IPDG")) {
    if(mesh->totalHaloPairs>0){
      ins->velocityHaloExtractKernel(mesh->Nelements,
                                 mesh->totalHaloPairs,
                                 mesh->o_haloElementList,
                                 ioffset,
                                 ins->o_U,
                                 ins->o_V,
                                 ins->o_vHaloBuffer);

      // copy extracted halo to HOST 
      ins->o_vHaloBuffer.copyTo(sendBuffer);           
    
      // start halo exchange
      meshHaloExchangeStart(mesh,
                           mesh->Np*(ins->NVfields)*sizeof(dfloat),
                           sendBuffer,
                           recvBuffer);
    }
  //}
  
  occaTimerTic(mesh->device,"DivergenceVolume");
  // computes div u^(n+1) volume term
  ins->divergenceVolumeKernel(mesh->Nelements,
                             mesh->o_vgeo,
                             mesh->o_DrT,
                             mesh->o_DsT,
                             ioffset,
                             ins->o_U,
                             ins->o_V,
                             ins->o_rhsP);
   occaTimerToc(mesh->device,"DivergenceVolume");

  //if (strstr(ins->pSolverOptions,"IPDG")) {
    if(mesh->totalHaloPairs>0){
      meshHaloExchangeFinish(mesh);

      ins->o_vHaloBuffer.copyFrom(recvBuffer); 

      ins->velocityHaloScatterKernel(mesh->Nelements,
                                    mesh->totalHaloPairs,
                                    mesh->o_haloElementList,
                                    ioffset,
                                    ins->o_U,
                                    ins->o_V,
                                    ins->o_vHaloBuffer);
    }

    occaTimerTic(mesh->device,"DivergenceSurface");
    //computes div u^(n+1) surface term
    ins->divergenceSurfaceKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_LIFTT,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                t,
                                mesh->o_x,
                                mesh->o_y,
                                ioffset,
                                ins->o_U,
                                ins->o_V,
                                ins->o_rhsP);
    occaTimerToc(mesh->device,"DivergenceSurface");
  //}

  
  occaTimerTic(mesh->device,"PoissonRhsForcing");
  // compute all forcing i.e. f^(n+1) - grad(Pr)
  ins->poissonRhsForcingKernel(mesh->Nelements,
                              mesh->o_vgeo,
                              mesh->o_MM,
                              ins->dt,  
                              ins->g0,
                              ins->o_rhsP);
  occaTimerToc(mesh->device,"PoissonRhsForcing");

#if 0
  //add penalty from jumps in previous pressure
  ins->poissonPenaltyKernel(mesh->Nelements,
                                mesh->o_sgeo,
                                mesh->o_vgeo,
                                mesh->o_DrT,
                                mesh->o_DsT,
                                mesh->o_LIFTT,
                                mesh->o_MM,
                                mesh->o_vmapM,
                                mesh->o_vmapP,
                                mesh->o_EToB,
                                ins->tau,
                                mesh->o_x,
                                mesh->o_y,
                                t,
                                ins->dt,
                                ins->c0,
                                ins->c1,
                                ins->c2,
                                ins->index,
                                (mesh->Nelements+mesh->totalHaloPairs),
                                ins->o_P,
                                ins->o_rhsP);
  #endif

  #if 1 // if time dependent BC
  //
  const int pressure_solve = 0; // ALGEBRAIC SPLITTING 
  if (strstr(ins->pSolverOptions,"CONTINUOUS")) {
    ins->poissonRhsBCKernel(mesh->Nelements,
                            pressure_solve,
                            mesh->o_ggeo,
                            mesh->o_sgeo,
                            mesh->o_SrrT,
                            mesh->o_SrsT,
                            mesh->o_SsrT,
                            mesh->o_SssT,
                            mesh->o_MM,
                            mesh->o_vmapM,
                            mesh->o_sMT,
                            t,
                            ins->dt,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_mapB,
                            ins->o_rhsP);
  } else if (strstr(ins->pSolverOptions,"IPDG")) {
    occaTimerTic(mesh->device,"PoissonRhsIpdg"); 
    ins->poissonRhsIpdgBCKernel(mesh->Nelements,
                                  pressure_solve,
                                  mesh->o_vmapM,
                                  mesh->o_vmapP,
                                  ins->tau,
                                  t,
                                  ins->dt,
                                  mesh->o_x,
                                  mesh->o_y,
                                  mesh->o_vgeo,
                                  mesh->o_sgeo,
                                  mesh->o_EToB,
                                  mesh->o_DrT,
                                  mesh->o_DsT,
                                  mesh->o_LIFTT,
                                  mesh->o_MM,
                                  ins->o_rhsP);
    occaTimerToc(mesh->device,"PoissonRhsIpdg");
  }
  #endif

  // gather-scatter
  if(strstr(ins->pSolverOptions, "CONTINUOUS")){
    ellipticParallelGatherScatterTri2D(mesh, mesh->ogs, ins->o_rhsP, ins->o_rhsP, dfloatString, "add");  
    if (mesh->Nmasked) mesh->maskKernel(mesh->Nmasked, mesh->o_maskIds, ins->o_rhsP);
  }


  if(!ins->maxPresHistory) { //if not storing successive pressure fields just solve

    occaTimerTic(mesh->device,"Pr Solve");
    ins->NiterP = ellipticSolveTri2D(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI,  ins->pSolverOptions); 
    occaTimerToc(mesh->device,"Pr Solve"); 
  
  } else {
    int izero = 0;
    dfloat zero = 0., one = 1.0;

    dfloat *blockReduction = ins->blockReduction;
    int Nblock = ins->Nblock;
    int Ntotal = mesh->Nelements*mesh->Np;
    occa::memory &o_blockReduction = ins->o_blockReduction;

    if (ins->NpresHistory) { //if we have at least one pressure history use projection
      // inner product rhs with orthonormalized set of previous pressures
      // alpha_l = rhsP . pres_l
      if(strstr(ins->pSolverOptions,"CONTINUOUS"))
        ins->multiWeightedInnerProductKernel(ins->NpresHistory, Nblock, Ntotal, solver->o_invDegree, ins->o_presHistory, ins->o_rhsP, o_blockReduction);
      else
        ins->multiInnerProductKernel(ins->NpresHistory, Nblock, Ntotal, ins->o_presHistory, ins->o_rhsP, o_blockReduction);
      o_blockReduction.copyTo(blockReduction,ins->NpresHistory*Nblock*sizeof(dfloat));

      for (int l=0;l<ins->NpresHistory;l++) {
        ins->presAlpha[l] = 0.;
        ins->presLocalAlpha[l] = 0.;
        for(int n=0;n<Nblock;++n){
          ins->presLocalAlpha[l] += blockReduction[n+l*Nblock];
        }  
      }
      MPI_Allreduce(ins->presLocalAlpha, ins->presAlpha, ins->NpresHistory, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      ins->o_presAlpha.copyFrom(ins->presAlpha);

      // PIbar = sum_l alpha_l pres_l
      ins->multiScaledAddKernel(ins->NpresHistory, Ntotal, ins->o_presAlpha, ins->o_presHistory, zero, ins->o_PIbar);

      // APIbar = A* PIbar
      ellipticOperator2D(solver, 0.0, ins->o_PIbar, ins->o_APIbar, ins->pSolverOptions);

      // subtract r = b - A*x
      ellipticScaledAdd(solver, -1.f, ins->o_APIbar, 1.f, ins->o_rhsP);

      //solve
      occaTimerTic(mesh->device,"Pr Solve");
      ins->NiterP = ellipticSolveTri2D(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI,  ins->pSolverOptions); 
      occaTimerToc(mesh->device,"Pr Solve"); 

      // PI += PIbar
      ellipticScaledAdd(solver, 1.f, ins->o_PIbar, 1.f, ins->o_PI);
    } else { //just solve

      //solve
      occaTimerTic(mesh->device,"Pr Solve");
      ins->NiterP = ellipticSolveTri2D(solver, 0.0, ins->presTOL, ins->o_rhsP, ins->o_PI,  ins->pSolverOptions); 
      occaTimerToc(mesh->device,"Pr Solve"); 

    }


    //update orthonormal pressure history
    if (ins->NpresHistory == ins->maxPresHistory) { //max history already stored
      
      //reset 
      ins->NpresHistory = 1; 

      // APIbar = A* PI
      ellipticOperator2D(solver, 0.0, ins->o_PI, ins->o_APIbar, ins->pSolverOptions);

      //compute ||PI||_A
      if(strstr(ins->pSolverOptions,"CONTINUOUS"))
        ins->multiWeightedInnerProductKernel(ins->NpresHistory, Nblock, Ntotal, solver->o_invDegree, ins->o_PI, ins->o_APIbar, o_blockReduction);
      else
        ins->multiInnerProductKernel(ins->NpresHistory, Nblock, Ntotal, ins->o_PI, ins->o_APIbar, o_blockReduction);
      o_blockReduction.copyTo(blockReduction,Nblock*sizeof(dfloat));

      ins->presAlpha[0] = 0.;
      ins->presLocalAlpha[0] = 0.;
      for(int n=0;n<Nblock;++n){
        ins->presLocalAlpha[0] += blockReduction[n];
      }  
      MPI_Allreduce(ins->presLocalAlpha, ins->presAlpha, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      dfloat norm = 1/sqrt(ins->presAlpha[0]);

      // pres_1 = PI/||PI||_A
      ins->scaledAddKernel(Ntotal, norm, izero, ins->o_PI, zero, izero, ins->o_presHistory);
    
    } else { //add this solution to the orthonormal history

      if(ins->NpresHistory) { //if we have at least one pressure history
        // inner product rhs with orthonormalized set of previous pressures
        // alpha_l = - PI . pres_l
        if(strstr(ins->pSolverOptions,"CONTINUOUS"))
          ins->multiWeightedInnerProductKernel(ins->NpresHistory, Nblock, Ntotal, solver->o_invDegree, ins->o_presHistory, ins->o_PI, o_blockReduction);
        else
          ins->multiInnerProductKernel(ins->NpresHistory, Nblock, Ntotal, ins->o_presHistory, ins->o_PI, o_blockReduction);
        o_blockReduction.copyTo(blockReduction,ins->NpresHistory*Nblock*sizeof(dfloat));

        for (int l=0;l<ins->NpresHistory;l++) {
          ins->presAlpha[l] = 0.;
          ins->presLocalAlpha[l] = 0.;
          for(int n=0;n<Nblock;++n){
            ins->presLocalAlpha[l] -= blockReduction[n+l*Nblock];
          }  
        }
        MPI_Allreduce(ins->presLocalAlpha, ins->presAlpha, ins->NpresHistory, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
        ins->o_presAlpha.copyFrom(ins->presAlpha);

        // PIbar = PI + sum_l alpha_l pres_l
        ins->o_PIbar.copyFrom(ins->o_PI,Ntotal*sizeof(dfloat));
        ins->multiScaledAddKernel(ins->NpresHistory, Ntotal, ins->o_presAlpha, ins->o_presHistory, one, ins->o_PIbar);
      } else {

        //just copy
        ins->o_PIbar.copyFrom(ins->o_PI,Ntotal*sizeof(dfloat));
      }

      // APIbar = A* PIbar
      ellipticOperator2D(solver, 0.0, ins->o_PIbar, ins->o_APIbar, ins->pSolverOptions);

      //compute ||PIbar||_A
      if(strstr(ins->pSolverOptions,"CONTINUOUS"))
        ins->multiWeightedInnerProductKernel(1, Nblock, Ntotal, solver->o_invDegree, ins->o_PIbar, ins->o_APIbar, o_blockReduction);
      else
        ins->multiInnerProductKernel(1, Nblock, Ntotal, ins->o_PI, ins->o_APIbar, o_blockReduction);
      o_blockReduction.copyTo(blockReduction,Nblock*sizeof(dfloat));

      ins->presAlpha[0] = 0.;
      ins->presLocalAlpha[0] = 0.;
      for(int n=0;n<Nblock;++n){
        ins->presLocalAlpha[0] += blockReduction[n];
      }  
      MPI_Allreduce(ins->presLocalAlpha, ins->presAlpha, 1, MPI_DFLOAT, MPI_SUM, MPI_COMM_WORLD);
      dfloat norm = 1/sqrt(ins->presAlpha[0]);
   
      // pres_1 = PI/||PI||_A
      ins->scaledAddKernel(Ntotal, norm, izero, ins->o_PIbar, zero, ins->NpresHistory*Ntotal, ins->o_presHistory);
    
      ins->NpresHistory++;
    }
  }

  if (strstr(ins->pSolverOptions,"CONTINUOUS")) {
    ins->poissonAddBCKernel(mesh->Nelements,
                            pressure_solve,
                            t,
                            ins->dt,
                            mesh->o_x,
                            mesh->o_y,
                            mesh->o_mapB,
                            ins->o_PI);
  }
}
