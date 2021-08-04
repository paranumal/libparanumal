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

#include "lbs.hpp"

// interpolate data to plot nodes and save to file (one per process)
void lbs_t::PlotFields(dfloat* Q, dfloat *V, char *fileName){

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

  //scratch space for interpolation
  size_t NscratchBytes = mymax(mesh.Np, mesh.plotNp)*sizeof(dfloat);
  dfloat* scratch = (dfloat *) malloc(2*NscratchBytes);

  dfloat* Ix = (dfloat *) malloc(mesh.plotNp*sizeof(dfloat));
  dfloat* Iy = (dfloat *) malloc(mesh.plotNp*sizeof(dfloat));
  dfloat* Iz = (dfloat *) malloc(mesh.plotNp*sizeof(dfloat));

  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh.Nelements;++e){
    mesh.PlotInterp(mesh.x + e*mesh.Np, Ix, scratch);
    mesh.PlotInterp(mesh.y + e*mesh.Np, Iy, scratch);
    mesh.PlotInterp(mesh.z + e*mesh.Np, Iz, scratch);

    for(int n=0;n<mesh.plotNp;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", Ix[n],Iy[n],Iz[n]);
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  free(Ix); free(Iy); free(Iz);

  dfloat* Ir = (dfloat *) malloc(mesh.plotNp*sizeof(dfloat));
  dfloat* Iu = (dfloat *) malloc(mesh.plotNp*sizeof(dfloat));
  dfloat* Iv = (dfloat *) malloc(mesh.plotNp*sizeof(dfloat));
  dfloat* Iw = (dfloat *) malloc(mesh.plotNp*sizeof(dfloat));

  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  if (U!=nullptr) {
    // write out velocity
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Velocity\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", mesh.dim);
    for(dlong e=0;e<mesh.Nelements;++e){
      mesh.PlotInterp(U + 1*mesh.Np + e*mesh.Np*Nmacro, Iu, scratch);
      mesh.PlotInterp(U + 2*mesh.Np + e*mesh.Np*Nmacro, Iv, scratch);
      if(mesh.dim==3)
        mesh.PlotInterp(U + 3*mesh.Np + e*mesh.Np*Nmacro, Iw, scratch);

      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        fprintf(fp, "       ");
        if (mesh.dim==2)
          fprintf(fp, "%g %g\n", Iu[n], Iv[n]);
        else
          fprintf(fp, "%g %g %g\n", Iu[n], Iv[n], Iw[n]);
      }
    }
    fprintf(fp, "       </DataArray>\n");
  }

  if (U!=nullptr) {
    // write out pressure
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Density\" Format=\"ascii\">\n");
    for(dlong e=0;e<mesh.Nelements;++e){
      mesh.PlotInterp(U + 0*mesh.Np + e*mesh.Np*Nmacro, Ir, scratch);

      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        fprintf(fp, "%g\n", Ir[n]);
      }
    }
    fprintf(fp, "       </DataArray>\n");
  }

  if (V!=nullptr) {
    // write out vorticity
    if(mesh.dim==2){
      fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" Format=\"ascii\">\n");
      for(dlong e=0;e<mesh.Nelements;++e){
        mesh.PlotInterp(V + e*mesh.Np, Ir, scratch);

        for(int n=0;n<mesh.plotNp;++n){
          fprintf(fp, "       ");
          fprintf(fp, "%g\n", Ir[n]);
        }
      }
    } else {
      fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Vorticity\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
      for(dlong e=0;e<mesh.Nelements;++e){
        mesh.PlotInterp(V + 0*mesh.Np + e*mesh.Np*3, Iu, scratch);
        mesh.PlotInterp(V + 1*mesh.Np + e*mesh.Np*3, Iv, scratch);
        mesh.PlotInterp(V + 2*mesh.Np + e*mesh.Np*3, Iw, scratch);

        for(int n=0;n<mesh.plotNp;++n){
          fprintf(fp, "       ");
          fprintf(fp, "       ");
          fprintf(fp, "%g %g %g\n", Iu[n], Iv[n], Iw[n]);
        }
      }
    }
    fprintf(fp, "       </DataArray>\n");
  }
  fprintf(fp, "     </PointData>\n");

  free(Ir); free(Iu); free(Iv); free(Iw);

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

  free(scratch);
}
