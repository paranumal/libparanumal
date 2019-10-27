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

#include "cns.hpp"

// interpolate data to plot nodes and save to file (one per process)
void cns_t::PlotFields(dfloat* Q, dfloat *V, char *fileName){

  FILE *fp;

  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          mesh.Nelements*mesh.plotNp,
          mesh.Nelements*mesh.plotNelements);

  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.plotNp;++n){
      dfloat plotxn = 0, plotyn = 0, plotzn = 0;

      for(int m=0;m<mesh.Np;++m){
        plotxn += mesh.plotInterp[n*mesh.Np+m]*mesh.x[m+e*mesh.Np];
        plotyn += mesh.plotInterp[n*mesh.Np+m]*mesh.y[m+e*mesh.Np];
        plotzn += mesh.plotInterp[n*mesh.Np+m]*mesh.z[m+e*mesh.Np];
      }

      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");


  if (Q!=NULL) {
    // write out density
    fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.plotNp;++n){
        dfloat plotpn = 0;
        for(int m=0;m<mesh.Np;++m){
          dfloat pm = Q[e*mesh.Np*Nfields+m];
          plotpn += mesh.plotInterp[n*mesh.Np+m]*pm;
        }

        fprintf(fp, "       ");
        fprintf(fp, "%g\n", plotpn);
      }
    }
    fprintf(fp, "       </DataArray>\n");

    // write out velocity
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", mesh.dim);
    for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.plotNp;++n){
        dfloat plotun = 0, plotvn = 0, plotwn = 0;
        for(int m=0;m<mesh.Np;++m){
          dfloat rm = Q[e*mesh.Np*Nfields+m+mesh.Np*0];
          dfloat um = Q[e*mesh.Np*Nfields+m+mesh.Np*1]/rm;
          dfloat vm = Q[e*mesh.Np*Nfields+m+mesh.Np*2]/rm;

          plotun += mesh.plotInterp[n*mesh.Np+m]*um;
          plotvn += mesh.plotInterp[n*mesh.Np+m]*vm;

          if(mesh.dim==3){
            dfloat wm = Q[e*mesh.Np*Nfields+m+mesh.Np*3]/rm;

            plotwn += mesh.plotInterp[n*mesh.Np+m]*wm;
          }
        }

        fprintf(fp, "       ");
        fprintf(fp, "       ");
        if (mesh.dim==2)
          fprintf(fp, "%g %g\n", plotun, plotvn);
        else
          fprintf(fp, "%g %g %g\n", plotun, plotvn, plotwn);
      }
    }
    fprintf(fp, "       </DataArray>\n");

    if (!isothermal) {
      const int eID = (mesh.dim==3) ? 4:3;
      // write out pressure
      fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Pressure\" Format=\"ascii\">\n");
      for(dlong e=0;e<mesh.Nelements;++e){
        for(int n=0;n<mesh.plotNp;++n){
          dfloat plotpn = 0;
          for(int m=0;m<mesh.Np;++m){
            dfloat rm = Q[e*mesh.Np*Nfields+m+mesh.Np*0];
            dfloat um = Q[e*mesh.Np*Nfields+m+mesh.Np*1]/rm;
            dfloat vm = Q[e*mesh.Np*Nfields+m+mesh.Np*2]/rm;
            dfloat wm = (mesh.dim==3) ? Q[e*mesh.Np*Nfields+m+mesh.Np*3]/rm : 0.0;
            dfloat em = Q[e*mesh.Np*Nfields+m+mesh.Np*eID];

            dfloat pm = (gamma-1)*(em-0.5*rm*(um*um+vm*vm+wm*wm));
            plotpn += mesh.plotInterp[n*mesh.Np+m]*pm;
          }

          fprintf(fp, "       ");
          fprintf(fp, "%g\n", plotpn);
        }
      }
      fprintf(fp, "       </DataArray>\n");
    }
  }

  if (V!=NULL) {
    // write out vorticity
    if(mesh.dim==2){
      fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" Format=\"ascii\">\n");
      for(dlong e=0;e<mesh.Nelements;++e){
        for(int n=0;n<mesh.plotNp;++n){
          dfloat plotVort = 0;
          for(int m=0;m<mesh.Np;++m){
            dlong id = m+e*mesh.Np;
            dfloat vort = V[id];
            plotVort += mesh.plotInterp[n*mesh.Np+m]*vort;
          }

          fprintf(fp, "       ");
          fprintf(fp, "%g\n", plotVort);
        }
      }
    } else {
      fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
      for(dlong e=0;e<mesh.Nelements;++e){
        for(int n=0;n<mesh.plotNp;++n){
          dfloat plotVortx = 0, plotVorty = 0, plotVortz = 0;
          for(int m=0;m<mesh.Np;++m){
            dlong id = m+e*mesh.Np*3;
            dfloat vortx = V[id+mesh.Np*0];
            dfloat vorty = V[id+mesh.Np*1];
            dfloat vortz = V[id+mesh.Np*2];
            plotVortx += mesh.plotInterp[n*mesh.Np+m]*vortx;
            plotVorty += mesh.plotInterp[n*mesh.Np+m]*vorty;
            plotVortz += mesh.plotInterp[n*mesh.Np+m]*vortz;
          }

          fprintf(fp, "       ");
          fprintf(fp, "%g %g %g\n", plotVortx, plotVorty, plotVortz);
        }
      }
    }
    fprintf(fp, "       </DataArray>\n");
  }

  fprintf(fp, "     </PointData>\n");

  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");

  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.plotNelements;++n){
      fprintf(fp, "       ");
      for(int m=0;m<mesh.plotNverts;++m){
        fprintf(fp, "%d ", e*mesh.plotNp + mesh.plotEToV[n*mesh.plotNverts+m]);
      }
      fprintf(fp, "\n");
    }
  }
  fprintf(fp, "        </DataArray>\n");

  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.plotNelements;++n){
      cnt += mesh.plotNverts;
      fprintf(fp, "       ");
      fprintf(fp, "%d\n", cnt);
    }
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.plotNelements;++n){
      if(mesh.dim==2)
        fprintf(fp, "5\n");
      else
        fprintf(fp, "10\n");
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

}
