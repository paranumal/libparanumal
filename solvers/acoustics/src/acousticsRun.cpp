/*

The MIT License (MIT)

Copyright (c) 2017-2022 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

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

#include "acoustics.hpp"
#include "timer.hpp"

void acoustics_t::Run(){

  dfloat startTime, finalTime;
  settings.getSetting("START TIME", startTime);
  settings.getSetting("FINAL TIME", finalTime);

  initialConditionKernel(mesh.Nelements,
                         startTime,
                         mesh.o_x,
                         mesh.o_y,
                         mesh.o_z,
                         o_q);

  dfloat cfl=1.0;
  settings.getSetting("CFL NUMBER", cfl);

  // set time step
  dfloat hmin = mesh.MinCharacteristicLength();
  dfloat vmax = MaxWaveSpeed();

  dfloat dt = cfl*hmin/(vmax*(mesh.N+1.)*(mesh.N+1.));
  timeStepper.SetTimeStep(dt);

  timePoint_t starts = GlobalPlatformTime(platform);

  timeStepper.Run(*this, o_q, startTime, finalTime);

  timePoint_t ends = GlobalPlatformTime(platform);

  double elapsedTime = ElapsedTime(starts, ends);
  std::cout << "elapsedTime = " << std::scientific << elapsedTime << std::endl;

  // output norm of final solution
  {
    //compute q.M*q
    dlong Nentries = mesh.Nelements*mesh.Np*Nfields;
    deviceMemory<dfloat> o_Mq = platform.reserve<dfloat>(Nentries);
    
    mesh.MassMatrixApply(o_q, o_Mq);

    dfloat norm2 = sqrt(platform.linAlg().innerProd(Nentries, o_q, o_Mq, mesh.comm));

    if(mesh.rank==0)
      printf("Solution norm = %17.15lg\n", norm2);

    deviceMemory<dfloat> o_exactqT = platform.reserve<dfloat>(Nentries);
    initialConditionKernel(mesh.Nelements,
                           finalTime,
                           mesh.o_x,
                           mesh.o_y,
                           mesh.o_z,
                           o_exactqT);
    
    memory<dfloat> exactqT(Nentries);
    o_exactqT.copyTo(exactqT);
    
    o_q.copyTo(q);

    dfloat errPLMax = 0;
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Np;++n){
        dlong id = e*mesh.Np*Nfields + n;
        dfloat errPLn = fabs(q[id] - exactqT[id]);
        
        errPLMax = std::max(errPLMax, errPLn);
      }
    }
    
    std::cout << std::setprecision (6);
    std::cout << std::scientific;
    std::cout << " errPLMax = " << errPLMax << std::endl;

    {
      memory<dfloat> pL(mesh.Np*mesh.Nelements);
      
      for(dlong e=0;e<mesh.Nelements;++e){
        for(int n=0;n<mesh.Np;++n){
          pL[e*mesh.Np + n] = q[e*mesh.Np*Nfields+n];
        }
      }

      
      std::string name;
      settings.getSetting("OUTPUT FILE NAME", name);
      char fname[BUFSIZ];
      // write to binary file
      sprintf(fname, "SOLN_DP_%04d.bin",  mesh.rank);
      FILE *fp = fopen(fname, "w");
      fprintf(fp, "# libParanumal binary format: blocks of Np*Nel doubles for (x,y,PL,PL,..). First line is: dim,Np,Nel,Nfields,sizeof(dfloat)\n");
      fprintf(fp, "%d %d %d %d %d (x,y,PL,PL)\n", mesh.dim, mesh.Np, mesh.Nelements, 4, sizeof(dfloat));
      dlong NpNel = mesh.Np*mesh.Nelements;
      fwrite(mesh.x.ptr(), NpNel, sizeof(dfloat), fp);
      fwrite(mesh.y.ptr(), NpNel, sizeof(dfloat), fp);
      fwrite(pL.ptr(), NpNel, sizeof(dfloat), fp);
      fwrite(pL.ptr(), NpNel, sizeof(dfloat), fp);
      fclose(fp);
    }
  }
}
