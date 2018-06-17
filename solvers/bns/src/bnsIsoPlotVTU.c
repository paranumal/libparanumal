#include "bns.h"

// interpolate data to plot nodes and save to file (one per process
void bnsIsoPlotVTU(bns_t *bns, int isoNtris, dfloat *isoq, char *fileName){

  mesh_t *mesh = bns->mesh;

  FILE *fp;
  
  fp = fopen(fileName, "w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
          isoNtris*3,
          isoNtris);
  
  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");
  
  // compute plot node coordinates on the fly
  for(dlong n=0;n<isoNtris*3;++n){
    int id = n*(mesh->dim+bns->isoNfields);
    dfloat plotxn = isoq[id+0];
    dfloat plotyn = isoq[id+1];
    dfloat plotzn = isoq[id+2];
    
    fprintf(fp, "       ");
    fprintf(fp, "%g %g %g\n", plotxn,plotyn,plotzn);
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");
  

  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  // write out velocity
  fprintf(fp, "        <DataArray type=\"Float32\" Name=\"Q\" NumberOfComponents=\"%d\" Format=\"ascii\">\n", bns->isoNfields);
  for(dlong n=0;n<3*isoNtris;++n){

    fprintf(fp, "       ");
    for(int fld=0;fld<bns->isoNfields;++fld){
      int id = n*(mesh->dim+bns->isoNfields) + mesh->dim + fld;
      dfloat plotqfld = isoq[id];

      fprintf(fp, "%g ", plotqfld);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "       </DataArray>\n");

  fprintf(fp, "     </PointData>\n");
  
  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(dlong e=0;e<isoNtris;++e){
    fprintf(fp, "       ");
    for(int m=0;m<3;++m){
      fprintf(fp, "%d ", e*3+m);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  dlong cnt = 0;
  for(dlong e=0;e<isoNtris;++e){
    cnt += 3;
    fprintf(fp, "       ");
    fprintf(fp, "%d\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(dlong e=0;e<isoNtris;++e){
    fprintf(fp, "5\n"); // tri
  }
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);

}
