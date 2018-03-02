#include "afterparser.h"

void send_vtu(char *filename, data *content) {

  FILE *fp;
  fp = fopen(filename,"w");

  fprintf(fp, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">\n");
  fprintf(fp, "  <UnstructuredGrid>\n");
  fprintf(fp, "    <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", 
	  content->points,
	  content->elements);

  // write out nodes
  fprintf(fp, "      <Points>\n");
  fprintf(fp, "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n");

  for (iint n = 0; n < content->points; n++) {
    fprintf(fp, "       ");
    fprintf(fp, "%g %g %g\n", content->x[n],content->y[n],content->z[n]);
  }
  
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Points>\n");

  fprintf(fp, "      <PointData Scalars=\"scalars\">\n");
  
  for (int i = 0; i < 10; ++i) {
    fprintf(fp, "        <DataArray type=\"Float32\" Name=\"q_%d\" Format=\"ascii\">\n",i+1);
    for(iint n=0;n<content->points;++n){
      fprintf(fp, "       ");
      fprintf(fp, "%g\n", content->q[n*10+i]);
    }
    fprintf(fp, "       </DataArray>\n");
  }

  fprintf(fp, "    </PointData>\n");

  fprintf(fp, "    <Cells>\n");
  fprintf(fp, "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">\n");
  
  for(iint e=0;e<content->elements;++e){
    fprintf(fp, "       ");
    fprintf(fp, "%d %d %d %d\n", content->connectivity[e*4], content->connectivity[e*4+1],content->connectivity[e*4+2],content->connectivity[e*4+3]);
  }
  
  fprintf(fp, "        </DataArray>\n");
  
  fprintf(fp, "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">\n");
  iint cnt = 0;
  for(iint e=0;e<content->elements;++e){
    cnt += 4;
    fprintf(fp, "       ");
    fprintf(fp, "%d\n", cnt);
  }
  fprintf(fp, "       </DataArray>\n");
  
  fprintf(fp, "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">\n");
  for(iint e=0;e<content->elements;++e){
    fprintf(fp, "9\n");
  }

  
  fprintf(fp, "        </DataArray>\n");
  fprintf(fp, "      </Cells>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </UnstructuredGrid>\n");
  fprintf(fp, "</VTKFile>\n");

  fclose(fp);
}
