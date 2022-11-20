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

#include "esdg.hpp"

void esdg_t::saveCheckpoint(memory<dfloat> &outq, dfloat time){

  static int frame=0;

  if (settings.compareSetting("CHECKPOINT SAVE","TRUE")) {

    // output field files                                                                                                                                     
    std::string name;
    settings.getSetting("CHECKPOINT SAVE NAME", name);

    char fname[BUFSIZ];
    sprintf(fname, "%s_%d_%06d.bin", name.c_str(), mesh.rank, frame++);

    // only set up for raw data 
    FILE *fp = fopen(fname, "w");
    fprintf(fp, "%17.15lf %%%% time\n", time);
    fwrite(outq.ptr(), mesh.Np*Nfields*mesh.Nelements, sizeof(dfloat), fp);
    fclose(fp);
  }
}

void esdg_t::loadCheckpoint(memory<dfloat> &inq, dfloat &time){

  std::string name;
  settings.getSetting("CHECKPOINT LOAD NAME", name);
  
  // only set up for raw data 
  FILE *fp = fopen(name.c_str(), "r");
  fscanf(fp, "%lf", &time);
  char buf[BUFSIZ];
  fgets(buf, BUFSIZ, fp);
  
  fread(inq.ptr(), mesh.Np*Nfields*mesh.Nelements, sizeof(dfloat), fp);
  fclose(fp);
  
  
}
