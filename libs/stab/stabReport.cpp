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

#include "stab.hpp"

namespace libp {

void stab_t::Test(){
  printf("testing stab..... \n");

// create a field to test
memory<dfloat> qtest; 
qtest.malloc((mesh.Nelements + mesh.totalHaloPairs)*mesh.Np);
deviceMemory<dfloat> o_qtest; 
deviceMemory<dfloat> o_testRhs; 
o_qtest   = platform.malloc<dfloat>((mesh.Nelements + mesh.totalHaloPairs)*mesh.Np);
o_testRhs = platform.malloc<dfloat>((mesh.Nelements + mesh.totalHaloPairs)*mesh.Np);


for(int e=0; e<mesh.Nelements; e++){
  for(int n=0; n<mesh.Np; n++){
    qtest[e*mesh.Np + n] = mesh.x[e*mesh.Np + n]<0 ? 1 : 0.0; 
  }
}
o_qtest.copyFrom(qtest);

const dfloat time = 0.0; 
const int tstep   = 0; 
Detect(o_qtest, o_testRhs, time); 

Apply(o_qtest, o_testRhs, time); 

Report(time, tstep); 
printf("testing stab.....done \n");
}



void stab_t::Report(dfloat time, int tstep){

  static int frame=0;

  dlong Ndetected  = GetElementNumber(o_elementList); 

  if(mesh.rank==0)
    printf("%5.2f (%d), %d (time, timestep, # of detected elements)\n", time, tstep, Ndetected);

  if (settings.compareSetting("STAB OUTPUT TO FILE","TRUE")) {	
	  // copy data back to host
  	if(elementList.length()!=0){
      o_elementList.copyTo(elementList);
    }

     // detector field
    if(o_qdetector.length()!=0){ 
      o_qdetector.copyTo(qdetector);
    }
    // Artificial Viscosity
    if(o_viscosityActivation.length()!=0){ 
      o_viscosityActivation.copyTo(viscosityActivation);
    }

    if(o_viscosity.length()!=0){ 
      o_viscosity.copyTo(viscosity);
    }
    if(o_vertexViscosity.length()!=0){ 
      o_vertexViscosity.copyTo(vertexViscosity);
    }

    // Limiter

    // if(o_qv.length()!=0){ 
    //   o_qv.copyTo(qv);
    // }

    // if(o_qc.length()!=0){ 
    //   o_qc.copyTo(qc);
    // }

    //  if(o_DX.length()!=0){ 
    //   o_DX.copyTo(DX);
    // }


    {
      char fname[BUFSIZ];
      sprintf(fname, "detector_%04d_%04d.vtu", mesh.rank, frame);
      PlotElements(elementList, std::string(fname));
  
      sprintf(fname, "stabFields_%04d_%04d.vtu", mesh.rank, frame);    
      PlotFields(qdetector, fname);
    }
    
    frame++; 
  }
  
  
}


dlong stab_t::GetElementNumber(deviceMemory<dlong>& o_list){
dlong Nelm = platform.linAlg().sum(mesh.Nelements*Ndfields, o_list, comm); 
return Nelm;
  
}



// interpolate data to plot nodes and save to file (one per processor)
void stab_t::PlotFields(memory<dfloat> Q, const std::string fileName){

FILE *fp;

  fp = fopen(fileName.c_str(), "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          mesh.Nelements*mesh.plotNp,
          mesh.Nelements*mesh.plotNelements);

  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

  //scratch space for interpolation
  size_t Nscratch = std::max(mesh.Np, mesh.plotNp);
  memory<dfloat> scratch(2*Nscratch);

  memory<dfloat> Ix(mesh.plotNp);
  memory<dfloat> Iy(mesh.plotNp);
  memory<dfloat> Iz(mesh.plotNp);

  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh.Nelements;++e){
    mesh.PlotInterp(mesh.x + e*mesh.Np, Ix, scratch);
    mesh.PlotInterp(mesh.y + e*mesh.Np, Iy, scratch);
    if(mesh.dim==3)
      mesh.PlotInterp(mesh.z + e*mesh.Np, Iz, scratch);

    if (mesh.dim==2) {
      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        fprintf(fp, "%g %g %g\n", Ix[n],Iy[n],0.0);
      }
    } else {
      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        fprintf(fp, "%g %g %g\n", Ix[n],Iy[n],Iz[n]);
      }
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  memory<dfloat> Iu(mesh.plotNp);
  memory<dfloat> Iv(mesh.plotNp);
  memory<dfloat> Iw(mesh.plotNp);
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
 
  if (Q.length()!=0) {
    // write out the field
     fprintf(fp, "        <DataArray type=\"Float32\" Name=\"DetectField\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", Ndfields);
    for(dlong e=0;e<mesh.Nelements;++e){
      mesh.PlotInterp(Q + 0*mesh.Np + e*mesh.Np*Ndfields, Iu, scratch);
       if(Ndfields>1)
        mesh.PlotInterp(Q + 1*mesh.Np  + e*mesh.Np*Ndfields, Iv, scratch);

      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        if(Ndfields==1)
          fprintf(fp, "%g \n", Iu[n]);
        else if(Ndfields==2)
          fprintf(fp, "%g %g\n", Iu[n], Iv[n]);
      }
    }
    fprintf(fp, "       </DataArray>\n");
  }

 if (viscosity.length()!=0) {
     fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Viscosity\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", Ndfields);
    for(dlong e=0;e<mesh.Nelements;++e){

      mesh.PlotInterp(viscosity + 0*mesh.Np  + e*mesh.Np*Ndfields, Iu, scratch);

      if(Ndfields>1)
        mesh.PlotInterp(viscosity + 1*mesh.Np  + e*mesh.Np*Ndfields, Iv, scratch);
      if(Ndfields>2)
        mesh.PlotInterp(viscosity + 2*mesh.Np  + e*mesh.Np*Ndfields, Iw, scratch);

      for(int n=0;n<mesh.plotNp;++n){
        fprintf(fp, "       ");
        fprintf(fp, "       ");
        if(Ndfields==1)
          fprintf(fp, "%g \n", Iu[n]);
        else if(Ndfields==2)
          fprintf(fp, "%g %g\n", Iu[n], Iv[n]);
        else if(Ndfields==3)
          fprintf(fp, "%g %g %g\n", Iu[n], Iv[n], Iw[n]);
      }
    }
    fprintf(fp, "       </DataArray>\n");
  }



// if(vertexViscosity.length()!=0){
//  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"vertexVisc\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", Ndfields);
//   for(dlong e=0;e<mesh.Nelements;++e){
//       for(int n=0;n<mesh.Nverts;++n){
//           fprintf(fp, "       ");
//            for(int fld=0; fld<Ndfields; fld++){
//              fprintf(fp, "%g ", vertexViscosity[e*mesh.Nverts*Ndfields + fld*mesh.Nverts + n]);
//           }
//           fprintf(fp, "\n");
//       }
//     }
//   fprintf(fp, "       </DataArray>\n");

//   }




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



// interpolate data to plot nodes and save to file (one per processor)
void stab_t::PlotElements(memory<dlong> ElementList, const std::string fileName){

  FILE *fp;
  fp = fopen(fileName.c_str(), "w");
  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n",
          mesh.Nelements*mesh.Nverts, mesh.Nelements);

  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

  // compute plot node coordinates on the fly
  for(dlong e=0;e<mesh.Nelements;++e){
    for(int n=0;n<mesh.Nverts;++n){
      const int vnode = mesh.vertexNodes[n]; 
      if(mesh.dim==2){
      fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", mesh.x[e*mesh.Np + vnode], 
                    mesh.y[e*mesh.Np + vnode],
                    0.0);
      }else{
        fprintf(fp, "       ");
      fprintf(fp, "%g %g %g\n", mesh.x[e*mesh.Np + vnode], 
                    mesh.y[e*mesh.Np + vnode],
                    mesh.z[e*mesh.Np + vnode]);
      }
    }
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");


  // write out field
  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  if(ElementList.length()!=0){
   fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Elements\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", Ndfields);
  for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Nverts;++n){
          fprintf(fp, "       ");
           for(int fld=0; fld<Ndfields; fld++){
             fprintf(fp, "%g ", dfloat(ElementList[e*Ndfields + fld]));
          }
          fprintf(fp, "\n");
      }
    }
  fprintf(fp, "       </DataArray>\n");
  }


if(viscosityActivation.length()!=0){
 fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Viscous Activation\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", Ndfields);
  for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Nverts;++n){
          fprintf(fp, "       ");
           for(int fld=0; fld<Ndfields; fld++){
             fprintf(fp, "%g ", viscosityActivation[e*Ndfields + fld]);
          }
          fprintf(fp, "\n");
      }
    }
  fprintf(fp, "       </DataArray>\n");
  }


  if(vertexViscosity.length()!=0){
 fprintf(fp, "        <DataArray type=\"Float32\" Name=\"vertexVisc\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", Ndfields);
  for(dlong e=0;e<mesh.Nelements;++e){
      for(int n=0;n<mesh.Nverts;++n){
          fprintf(fp, "       ");
           for(int fld=0; fld<Ndfields; fld++){
             fprintf(fp, "%g ", vertexViscosity[e*mesh.Nverts*Ndfields + fld*mesh.Nverts + n]);
          }
          fprintf(fp, "\n");
      }
    }
  fprintf(fp, "       </DataArray>\n");

  }








  fprintf(fp, "     </PointData>\n");
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh.Nelements;++e){
    fprintf(fp, "       ");
    for(int m=0;m<mesh.Nverts;++m){
      fprintf(fp, "%d ", e*mesh.Nverts + m);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "        </DataArray>\n");

  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong e=0;e<mesh.Nelements;++e){
    cnt += mesh.Nverts;
    fprintf(fp, "       ");
    fprintf(fp, "%d\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");

  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<mesh.Nelements;++e){
    if(mesh.elementType==Mesh::TRIANGLES)
      fprintf(fp, "5\n");
    else if(mesh.elementType==Mesh::TETRAHEDRA)
      fprintf(fp, "10\n");
    else if(mesh.elementType==Mesh::QUADRILATERALS)
      fprintf(fp, "9\n");
    else if(mesh.elementType==Mesh::HEXAHEDRA)
      fprintf(fp, "12\n");

  }
  
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
  
}








}