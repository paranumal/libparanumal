/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

#include "ins.h"
void insRestartWrite(ins_t *ins, setupAide &options, dfloat t){

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
  sprintf(fname, "%s_%04d.dat",(char*)outName.c_str(), mesh->rank);
  FILE *fp;
  
  // Overwrite  Binary File
  fp = fopen(fname,"wb");
  
  // Write q1....qN and all pml variables
  dfloat *elmField  = (dfloat *) calloc((ins->NVfields+1), sizeof(dfloat)); 
  dfloat *elmField2 = (dfloat *) calloc(2*ins->NVfields, sizeof(dfloat)); 

  // First Write the Solution Time
  fwrite(&t, sizeof(dfloat), 1, fp);
  //Write dt
  fwrite(&ins->dt, sizeof(dfloat), 1, fp);
  // Write output frame to prevent overwriting vtu files
  fwrite(&ins->frame, sizeof(int), 1, fp);

  // 
 if(options.compareArgs("TIME INTEGRATOR", "EXTBDF") ){

  // Write U and P 
  for(int s =0; s<ins->Nstages; s++){
    for(dlong e = 0;e<mesh->Nelements; e++){
      for(int n=0; n<mesh->Np; n++ ){
        const dlong idv = e*mesh->Np + n + s*ins->fieldOffset*ins->NVfields; 
        const dlong idp = e*mesh->Np + n + s*ins->fieldOffset; 
          for(int vf = 0; vf<ins->NVfields; vf++){
            elmField[vf]   =  ins->U[idv + vf*ins->fieldOffset];
          }
          elmField[ins->NVfields] = ins->P[idp]; 
        fwrite(elmField, sizeof(dfloat), (ins->NVfields+1), fp);
      }
    } 
  }
  // Write nonlinear History 
  for(int s =0; s<ins->Nstages; s++){
    for(dlong e = 0;e<mesh->Nelements; e++){
      for(int n=0; n<mesh->Np; n++ ){
        const dlong idv = e*mesh->Np + n + s*ins->fieldOffset*ins->NVfields; 
          
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

  mesh_t *mesh = ins->mesh;
  
  // Create Binary File Name
  char fname[BUFSIZ];
  string outName;
  options.getArgs("RESTART FILE NAME", outName);
  sprintf(fname, "%s_%04d.dat",(char*)outName.c_str(), mesh->rank);
  FILE *fp; 


  ins->restartedFromFile = 0; 
  // Overwrite  Binary File
  fp = fopen(fname,"rb");

  if(fp != NULL){
  
      // Write q1....qN and all pml variables
    dfloat *elmField  = (dfloat *) calloc((ins->NVfields+1), sizeof(dfloat)); 
    dfloat *elmField2 = (dfloat *) calloc(2*ins->NVfields, sizeof(dfloat)); 

    // First Write the Solution Time
    dfloat startTime = 0.0, dtold = 0.0;  
    fread(&startTime, sizeof(dfloat), 1, fp); 
    //Write dt
    fread(&dtold, sizeof(dfloat), 1, fp);
    // Write output frame to prevent overwriting vtu files
    fread(&ins->frame, sizeof(int), 1, fp);

    // 
    if(options.compareArgs("TIME INTEGRATOR", "EXTBDF") ){

      // Write U and P
      for(int s =0; s<ins->Nstages; s++){
        for(dlong e = 0;e<mesh->Nelements; e++){
          for(int n=0; n<mesh->Np; n++ ){

            fread(elmField, sizeof(dfloat), (ins->NVfields+1), fp);

            const dlong idv = e*mesh->Np + n + s*ins->fieldOffset*ins->NVfields; 
            const dlong idp = e*mesh->Np + n + s*ins->fieldOffset;

              for(int vf = 0; vf<ins->NVfields; vf++)
                ins->U[idv + vf*ins->fieldOffset] = elmField[vf];

              ins->P[idp] = elmField[ins->NVfields]; 
          }
        } 
      }
      // Write nonlinear History 
      for(int s =0; s<ins->Nstages; s++){
        for(dlong e = 0;e<mesh->Nelements; e++){
          for(int n=0; n<mesh->Np; n++ ){
            const dlong idv = e*mesh->Np + n + s*ins->fieldOffset*ins->NVfields;       
            fread(elmField2, sizeof(dfloat), 2*ins->NVfields, fp);       
              for(int vf = 0; vf<ins->NVfields; vf++)
                ins->NU[idv + vf*ins->fieldOffset] = elmField2[vf];    
              for(int vf = 0; vf<ins->NVfields; vf++)
                ins->GP[idv + vf*ins->fieldOffset] = elmField2[ins->NVfields+vf];
          }
        } 
      }
    }else{

      if(mesh->rank==0) printf("restart for ARK has not tested yet\n");
    }

  fclose(fp);

  ins->restartedFromFile = 1;  
  // Just Update start time
  ins->startTime = startTime; 
  ins->dt        = ins->dti; // set time-step to initial time-step estimate

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

   // Interpolate history if dt is changed in new setup
    if(fabs(ins->dt - dtold) > 1.E-12){ // guard for infinitesmall dt change

      if(mesh->rank==0) printf("Interpolating history for new dt...");

      // ins->Nstages is fixed to max 3; 
      dfloat interp[ins->Nstages][ins->Nstages]; 

      const dfloat tn0 = -0*dtold;
      const dfloat tn1 = -1*dtold;
      const dfloat tn2 = -2*dtold;

      for(int stage = 0; stage<ins->Nstages; stage++){

        dfloat t = -stage*ins->dt; // current time levels to be interpolated

        switch(ins->Nstages){
          case 1:
            interp[stage][0] = 1.f; interp[stage][1] = 0.f; interp[stage][2] = 0.f;
          break;
          case 2:
            interp[stage][0] = (t-tn1)/(tn0-tn1);
            interp[stage][1] = (t-tn0)/(tn1-tn0);
            interp[stage][2] = 0.f; 
          break;
          case 3:
            interp[stage][0] = (t-tn1)*(t-tn2)/((tn0-tn1)*(tn0-tn2)); 
            interp[stage][1] = (t-tn0)*(t-tn2)/((tn1-tn0)*(tn1-tn2));
            interp[stage][2] = (t-tn0)*(t-tn1)/((tn2-tn0)*(tn2-tn1));
          break;
        }
      }

      // Interpolation Storage after stage update 
      dfloat *NUs = (dfloat *)calloc(ins->NVfields*ins->Nstages,sizeof(dfloat));
      dfloat *GPs = (dfloat *)calloc(ins->NVfields*ins->Nstages,sizeof(dfloat));
      dfloat *Us  = (dfloat *)calloc(ins->NVfields*ins->Nstages,sizeof(dfloat));
      dfloat *Ps  = (dfloat *)calloc(ins->Nstages,sizeof(dfloat));

      // hold for all fileds
      dfloat *ui   = (dfloat *)calloc(ins->NVfields,sizeof(dfloat));      
      dfloat *nui  = (dfloat *)calloc(ins->NVfields,sizeof(dfloat));      
      dfloat *gpi  = (dfloat *)calloc(ins->NVfields,sizeof(dfloat));      
      
      // For EXTBDF history fields that have to be updated, NU, GP
      for(dlong e=0; e<mesh->Nelements; e++){
        for(int n=0; n<mesh->Np; n++){
         
          const dlong id           = e*mesh->Np + n; 
          const dlong stageOffset  = ins->NVfields*ins->fieldOffset;

          for(int stage = 0; stage<ins->Nstages; stage++){
            
            // Initialize to zero
            for(int fld = 0; fld<ins->NVfields; fld++){
              ui[fld]  = 0.0; nui[fld] = 0.0;  gpi[fld] = 0.0; 
            }

            dfloat pi = 0.0; 

            for(int s = 0; s<ins->Nstages; s++){
              for(int fld = 0; fld<ins->NVfields; fld++){
                 ui[fld]  += interp[stage][s]*ins->U [id+fld*ins->fieldOffset+s*stageOffset];
                nui[fld]  += interp[stage][s]*ins->NU[id+fld*ins->fieldOffset+s*stageOffset];
                gpi[fld]  += interp[stage][s]*ins->GP[id+fld*ins->fieldOffset+s*stageOffset];
              }
              pi += interp[stage][s]*ins->P[id+s*ins->fieldOffset];
            }

            // update field for this stage i.e. t = -stage*ins->dtNew
            for(int fld = 0; fld<ins->NVfields; fld++){
               Us[stage*ins->NVfields + fld] =  ui[fld]; 
              NUs[stage*ins->NVfields + fld] = nui[fld]; 
              GPs[stage*ins->NVfields + fld] = gpi[fld]; 
            }
            Ps[stage] = pi;  
          }

          // Update node values 
          for(int stage = 0; stage<ins->Nstages; stage++){
            for(int fld = 0; fld<ins->NVfields; fld++){
             ins-> U[id+fld*ins->fieldOffset+stage*stageOffset] =   Us[stage*ins->NVfields + fld]; 
             ins->NU[id+fld*ins->fieldOffset+stage*stageOffset] =  NUs[stage*ins->NVfields + fld]; 
             ins->GP[id+fld*ins->fieldOffset+stage*stageOffset] =  GPs[stage*ins->NVfields + fld]; 
             }
           ins->P[id + stage*ins->fieldOffset] = Ps[stage]; 
          }
        }
      }
    }

  }else{
      if(mesh->rank==0) printf("restart for ARK has not tested yet\n");
  }




  ins->o_U.copyFrom(ins->U);
  ins->o_P.copyFrom(ins->P);
  // History of nonlinear terms !!!!
  ins->o_NU.copyFrom(ins->NU);
  // History of Pressure !!!!!
  ins->o_GP.copyFrom(ins->GP);


}else{
  printf("No restart file...");
}

}





