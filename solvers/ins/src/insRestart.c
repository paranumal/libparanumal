#include "ins.h"
void insRestartWrite(ins_t *ins, setupAide &options, dfloat t){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = ins->mesh; 
  // Copy Field To Host
  // copy data back to host
  ins->o_U.copyTo(ins->U);
  ins->o_P.copyTo(ins->P);
  // History of nonlinear terms !!!!
  ins->o_NU.copyTo(ins->NU);
  // History of Pressure !!!!!
  ins->o_GP.copyTo(ins->GP);

  // Create Binary File Name
  char fname[BUFSIZ];
  string outName;
  options.getArgs("RESTART FILE NAME", outName);
  sprintf(fname, "%s_%04d.dat",(char*)outName.c_str(), rank);
  FILE *fp;
  
  // Overwrite  Binary File
  fp = fopen(fname,"wb");
  
  // Write q1....qN and all pml variables
  dfloat *elmField  = (dfloat *) calloc((ins->NVfields+1), sizeof(dfloat)); 
  dfloat *elmField2 = (dfloat *) calloc(2*ins->NVfields, sizeof(dfloat)); 

  // First Write the Solution Time
  fwrite(&t, sizeof(dfloat), 1, fp);
  // Write output frame to prevent overwriting vtu files
  fwrite(&ins->frame, sizeof(int), 1, fp);

  // 
 if(options.compareArgs("TIME INTEGRATOR", "EXTBDF") ){

  // Write U and P
  for(int s =0; s<ins->Nstages; s++){
    for(dlong e = 0;e<mesh->Nelements; e++){
      for(int n=0; n<mesh->Np; n++ ){
        const dlong idv = e*mesh->Np + n + s*ins->Ntotal*ins->NVfields; 
        const dlong idp = e*mesh->Np + n + s*ins->Ntotal; 
          for(int vf = 0; vf<ins->NVfields; vf++){
            elmField[vf]   =  ins->U[idv + vf*ins->fieldOffset];
          }
          elmField[ins->NVfields] = ins->P[idp]; 

        fwrite(elmField, sizeof(dfloat), (ins->NVfields+1), fp);
      }
    } 
  }
  // Write nonlinear History 
  for(int s =0; s<(ins->Nstages+1); s++){
    for(dlong e = 0;e<mesh->Nelements; e++){
      for(int n=0; n<mesh->Np; n++ ){
        const dlong idv = e*mesh->Np + n + s*ins->Ntotal*ins->NVfields; 
          
          for(int vf = 0; vf<ins->NVfields; vf++)
            elmField2[vf]   =  ins->NU[idv + vf*ins->fieldOffset];
          
          for(int vf = 0; vf<ins->NVfields; vf++)
            elmField2[ins->NVfields+ vf]   =  ins->GP[idv + vf*ins->fieldOffset];
        fwrite(elmField2, sizeof(dfloat), 2*ins->NVfields, fp);
      }
    } 
  }


}




fclose(fp); 


}




void insRestartRead(ins_t *ins, setupAide &options){

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  mesh_t *mesh = ins->mesh; 
  // Create Binary File Name
  char fname[BUFSIZ];
  string outName;
  options.getArgs("RESTART FILE NAME", outName);
  sprintf(fname, "%s_%04d.dat",(char*)outName.c_str(), rank);
  FILE *fp; 

  // Overwrite  Binary File
  fp = fopen(fname,"rb");

  if(fp != NULL){
  
      // Write q1....qN and all pml variables
    dfloat *elmField  = (dfloat *) calloc((ins->NVfields+1), sizeof(dfloat)); 
    dfloat *elmField2 = (dfloat *) calloc(2*ins->NVfields, sizeof(dfloat)); 

    // First Write the Solution Time
    dfloat startTime = 0.0; 
    fread(&startTime, sizeof(dfloat), 1, fp); 
    // Write output frame to prevent overwriting vtu files
    fread(&ins->frame, sizeof(int), 1, fp);

    // 
    if(options.compareArgs("TIME INTEGRATOR", "EXTBDF") ){

      // Write U and P
      for(int s =0; s<ins->Nstages; s++){
        for(dlong e = 0;e<mesh->Nelements; e++){
          for(int n=0; n<mesh->Np; n++ ){

            fread(elmField, sizeof(dfloat), (ins->NVfields+1), fp);

            const dlong idv = e*mesh->Np + n + s*ins->Ntotal*ins->NVfields; 
            const dlong idp = e*mesh->Np + n + s*ins->Ntotal; 
              for(int vf = 0; vf<ins->NVfields; vf++){
                ins->U[idv + vf*ins->fieldOffset] = elmField[vf];
              }
              ins->P[idp] = elmField[ins->NVfields]; 
          }
        } 
      }
      // Write nonlinear History 
      for(int s =0; s<(ins->Nstages+1); s++){
        for(dlong e = 0;e<mesh->Nelements; e++){
          for(int n=0; n<mesh->Np; n++ ){
            const dlong idv = e*mesh->Np + n + s*ins->Ntotal*ins->NVfields; 
            
            fread(elmField2, sizeof(dfloat), 2*ins->NVfields, fp);
              
              for(int vf = 0; vf<ins->NVfields; vf++)
                ins->NU[idv + vf*ins->fieldOffset] = elmField2[vf];
              
              for(int vf = 0; vf<ins->NVfields; vf++)
                ins->GP[idv + vf*ins->fieldOffset] = elmField2[ins->NVfields+vf];

          }
        } 
      }
    }else{

      printf("restart for ARK has not tested yet\n");


    }


  fclose(fp);

  ins->o_U.copyFrom(ins->U);
  ins->o_P.copyFrom(ins->P);
  // History of nonlinear terms !!!!
  ins->o_NU.copyFrom(ins->NU);
  // History of Pressure !!!!!
  ins->o_GP.copyFrom(ins->GP);

  // Just Update start time
  ins->startTime = startTime; 

  ins->dt = ins->dti; // set time-step to initial time-step estimate

   if (options.compareArgs("TIME INTEGRATOR", "EXTBDF")) {
    ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;

    if(ins->Nsubsteps){
      ins->dt         = ins->Nsubsteps*ins->dt;
      ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;
      ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
      ins->sdt        = ins->dt/ins->Nsubsteps;
    } else{
      ins->NtimeSteps = (ins->finalTime-ins->startTime)/ins->dt;
      ins->dt         = (ins->finalTime-ins->startTime)/ins->NtimeSteps;
    }
  }else{
      printf("restart for ARK has not tested yet\n");
  }



}

}