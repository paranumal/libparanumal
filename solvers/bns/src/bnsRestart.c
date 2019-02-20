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

#include "bns.h"
void bnsRestartWrite(bns_t *bns, setupAide &options, dfloat time){

  mesh_t *mesh = bns->mesh;
  // Copy Field To Host
  bns->o_q.copyTo(bns->q);

  // Create Binary File Name
  char fname[BUFSIZ];
  string outName;
  options.getArgs("RESTART FILE NAME", outName);
  sprintf(fname, "%s_%04d.dat",(char*)outName.c_str(), mesh->rank);
  FILE *fp;

  // Overwrite  Binary File
  fp = fopen(fname,"wb");

  // Write q1....qN and all pml variables
  dfloat *elmField  = (dfloat *) calloc(bns->Nfields, sizeof(dfloat));
  dfloat *pmlField  = (dfloat *) calloc(bns->Nfields*bns->dim, sizeof(dfloat));

  // First Write the Solution Time
  fwrite(&time, sizeof(dfloat), 1, fp);
  // Write output frame to prevent overwriting vtu files
  fwrite(&bns->frame, sizeof(int), 1, fp);

  //
 if(options.compareArgs("TIME INTEGRATOR", "LSERK") ||
    options.compareArgs("TIME INTEGRATOR", "SARK")){
  // Write only q works check for MRAB, write history
  for(dlong e = 0; e<mesh->Nelements; e++){
    for(int n=0; n<mesh->Np; n++ ){
      const dlong id = e*bns->Nfields*mesh->Np + n;
      for(int fld=0; fld<bns->Nfields; fld++){
        elmField[fld] =  bns->q[id + fld*mesh->Np];
      }
      fwrite(elmField, sizeof(dfloat), bns->Nfields, fp);
    }
  }

  if(mesh->pmlNelements){
    // Copy Field To Host
    bns->o_pmlqx.copyTo(bns->pmlqx);
    bns->o_pmlqy.copyTo(bns->pmlqy);
    if(bns->dim==3)
      bns->o_pmlqz.copyTo(bns->pmlqz);

    // Write node data
    for(dlong es = 0; es<mesh->pmlNelements; es++)
    {
      dlong e      = mesh->pmlElementIds[es];
      dlong pmlId  = mesh->pmlIds[es];
      for(int n=0; n<mesh->Np; n++ ){
        const dlong id = pmlId*bns->Nfields*mesh->Np + n;
        for(int fld=0; fld<bns->Nfields; fld++){
          pmlField[fld + 0*bns->Nfields] =  bns->pmlqx[id + fld*mesh->Np];
          pmlField[fld + 1*bns->Nfields] =  bns->pmlqy[id + fld*mesh->Np];
          if(bns->dim==3)
            pmlField[fld + 2*bns->Nfields] =  bns->pmlqz[id + fld*mesh->Np];
        }
      fwrite(pmlField, sizeof(dfloat), bns->Nfields*bns->dim, fp);
      }
    }
  }
}

#if 0
  // Assumes the same dt, there is no interpolation of history
  // will be mofied later if the time step size changes
   if(options.compareArgs("TIME INTEGRATOR", "MRSAAB")){
   // Keep the same index to prvent order loss at restart
    for(int l = 0; l<mesh->MRABNlevels; l++)
      fwrite(&mesh->MRABshiftIndex[l], sizeof(int),1, fp);

      // Copy right hand side history
      bns->o_rhsq.copyTo(bns->rhsq);

      const dlong offset   = mesh->Nelements*mesh->Np*bns->Nfields;

      // Write all history, order is not important
      // as long as reading with the same order
      for(int nrhs = 0; nrhs<bns->Nrhs; nrhs++){
        for(dlong e = 0; e<mesh->Nelements; e++){
          for(int n=0; n<mesh->Np; n++ ){
            const dlong id = e*bns->Nfields*mesh->Np + n + nrhs*offset;
            for(int fld=0; fld<bns->Nfields; fld++){
              elmField[fld] =  bns->q[id + fld*mesh->Np];
            }
          fwrite(elmField, sizeof(dfloat), bns->Nfields, fp);
          }
        }
      }

    if(bns->pmlFlag){
      // Copy Field To Host
      bns->o_pmlrhsqx.copyTo(bns->pmlrhsqx);
      bns->o_pmlrhsqy.copyTo(bns->pmlrhsqy);
      if(bns->dim==3)
        bns->o_pmlrhsqz.copyTo(bns->pmlrhsqz);

      const dlong pmloffset = mesh->pmlNelements*mesh->Np*bns->Nfields;

      // Write node data
      for(int nrhs = 0; nrhs<bns->Nrhs; nrhs++){
        for(dlong es = 0; es<mesh->pmlNelements; es++){
          dlong e      = mesh->pmlElementIds[es];
          dlong pmlId  = mesh->pmlIds[es];
          for(int n=0; n<mesh->Np; n++ ){
            const dlong id = pmlId*bns->Nfields*mesh->Np + n + nrhs*pmloffset;
            for(int fld=0; fld<bns->Nfields; fld++){
              pmlField[fld + 0*bns->Nfields] =  bns->pmlrhsqx[id + fld*mesh->Np];
              pmlField[fld + 1*bns->Nfields] =  bns->pmlrhsqy[id + fld*mesh->Np];
              if(bns->dim==3)
                pmlField[fld + 2*bns->Nfields] =  bns->pmlrhsqz[id + fld*mesh->Np];
              }
            fwrite(pmlField, sizeof(dfloat), bns->Nfields*bns->dim, fp);
          }
        }
      }
    }
   }
#endif

  fclose(fp);
}




void bnsRestartRead(bns_t *bns, setupAide &options){

  mesh_t *mesh = bns->mesh;
  // Create Binary File Name
  char fname[BUFSIZ];
  string outName;
  options.getArgs("RESTART FILE NAME", outName);
  sprintf(fname, "%s_%04d.dat",(char*)outName.c_str(), mesh->rank);
  FILE *fp;

  size_t nentries;

  // Overwrite  Binary File
  fp = fopen(fname,"rb");

  if(fp != NULL){

    // Write q1....qN and all pml variables
    dfloat *elmField  = (dfloat *) calloc(bns->Nfields, sizeof(dfloat));
    dfloat *pmlField  = (dfloat *) calloc(bns->Nfields*bns->dim, sizeof(dfloat));

    dfloat startTime = 0.0;
    // Update Start Time
    nentries = fread(&startTime, sizeof(dfloat),1, fp);
    // Update frame number to contioune outputs
    nentries = fread(&bns->frame, sizeof(int), 1, fp);


    if(mesh->rank==0) printf("Restart time: %.4e ...", startTime);

    // Write only q works check for MRAB, write history
    for(dlong e = 0; e<mesh->Nelements; e++){
      for(int n=0; n<mesh->Np; n++ ){

        const dlong id = e*bns->Nfields*mesh->Np + n;
        nentries = fread(elmField, sizeof(dfloat), bns->Nfields, fp);

        for(int fld=0; fld<bns->Nfields; fld++){
          bns->q[id + fld*mesh->Np] = elmField[fld];
        }

      }
    }

    if(bns->pmlFlag){

      // Write node data
      for(dlong es = 0; es<mesh->pmlNelements; es++)
      {
        dlong e      = mesh->pmlElementIds[es];
        dlong pmlId  = mesh->pmlIds[es];
        for(int n=0; n<mesh->Np; n++ ){
          const dlong id = pmlId*bns->Nfields*mesh->Np + n;

          nentries = fread(pmlField, sizeof(dfloat), bns->Nfields*bns->dim, fp);

          for(int fld=0; fld<bns->Nfields; fld++){
            bns->pmlqx[id + fld*mesh->Np] = pmlField[fld + 0*bns->Nfields];
            bns->pmlqy[id + fld*mesh->Np] = pmlField[fld + 1*bns->Nfields];
            if(bns->dim==3)
              bns->pmlqz[id + fld*mesh->Np] = pmlField[fld + 2*bns->Nfields];
          }
        }
      }
    }

  fclose(fp);

  // Update Fields
    bns->o_q.copyFrom(bns->q);

    if(mesh->pmlNelements){
      bns->o_pmlqx.copyFrom(bns->pmlqx);
      bns->o_pmlqy.copyFrom(bns->pmlqy);
      if(bns->dim==3)
        bns->o_pmlqz.copyFrom(bns->pmlqz);
    }

  // Just Update Time Step Size
  bns->startTime = startTime;
  // Update NtimeSteps and dt
   bns->NtimeSteps = (bns->finalTime-bns->startTime)/bns->dt;
   bns->dt         = (bns->finalTime-bns->startTime)/bns->NtimeSteps;
}else{
  printf("No restart file...");

}
}
