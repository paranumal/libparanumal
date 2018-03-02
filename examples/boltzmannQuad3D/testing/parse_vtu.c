#include "afterparser.h"

//filename is null-terminated
void parse_vtu(char *filename,data *parsed) {
  FILE *fp = fopen(filename, "r");
  
  char buf[BUFSIZ];

  //find xml line with point/element numbers
  fgets(buf,BUFSIZ,fp);
  while(!strstr(buf,"NumberOfPoints")) fgets(buf,BUFSIZ,fp);

  //extract point/element numbers
  char *buf_small = strtok(buf,"\"");
  parsed->points = atoi(strtok(NULL,"\""));
  buf_small = strtok(NULL,"\"");
  parsed->elements = atoi(strtok(NULL,"\""));

  //skip to header with point coordinates
  while(!strstr(buf,"<Points>")) fgets(buf,BUFSIZ,fp);

  //discard next header
  fgets(buf,BUFSIZ,fp);

  //load coordinates into array
  parsed->x = (dfloat *) calloc(parsed->points,sizeof(dfloat));
  parsed->y = (dfloat *) calloc(parsed->points,sizeof(dfloat));
  parsed->z = (dfloat *) calloc(parsed->points,sizeof(dfloat));

  for(int i = 0; i < parsed->points; ++i) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%lf %lf %lf",parsed->x+i,parsed->y+i,parsed->z+i);
  }

  //skip to start of mrab levels
  while(!strstr(buf,"mrab_levels")) fgets(buf,BUFSIZ,fp);

  parsed->mrab_levels = (iint *) calloc(parsed->elements,sizeof(iint));
  for(int i = 0; i < parsed->elements; i++) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%lf",parsed->mrab_levels+i);
    //all points are the same, so skip to next element
    for(int j = 0; j < parsed->points/parsed->elements - 1;j++)
      fgets(buf,BUFSIZ,fp);
  }

  //load all q values in a loop
  parsed->q = (dfloat *) calloc(parsed->points*10,sizeof(dfloat));
  
  for(int fld = 0; fld < 10; ++fld) {
    char cmp[5];
    sprintf(cmp,"q_%d",fld+1);
    while(!strstr(buf,cmp)) fgets(buf,BUFSIZ,fp);  
        
    for(int i = 0; i < parsed->points; i++) {
      fgets(buf,BUFSIZ,fp);
      sscanf(buf,"%lf",parsed->q+i*10+fld);
    }
  }

  //skip to start of density field
  while(!strstr(buf,"Density")) fgets(buf,BUFSIZ,fp);  

  parsed->density = (dfloat *) calloc(parsed->points,sizeof(dfloat));

  for(int i = 0; i < parsed->points; i++) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%lf",parsed->density+i);
  }

  //skip to start of jacobian
  while(!strstr(buf,"jacobian")) fgets(buf,BUFSIZ,fp);

  parsed->jacobian = (dfloat *) calloc(parsed->points,sizeof(dfloat));

  for(int i = 0; i < parsed->points; i++) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%lf",parsed->density+i);
  }

  while(!strstr(buf,"Velocity")) fgets(buf,BUFSIZ,fp);

  parsed->vel_x = (dfloat *) calloc(parsed->points,sizeof(dfloat));
  parsed->vel_y = (dfloat *) calloc(parsed->points,sizeof(dfloat));
  parsed->vel_z = (dfloat *) calloc(parsed->points,sizeof(dfloat));

  for(int i = 0; i < parsed->points; i++) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%lf %lf %lf",parsed->vel_x + i,parsed->vel_y+i,parsed->vel_z+i);
  }

  while(!strstr(buf,"Vorticity")) fgets(buf,BUFSIZ,fp);

  parsed->vort_x = (dfloat *) calloc(parsed->points,sizeof(dfloat));
  parsed->vort_y = (dfloat *) calloc(parsed->points,sizeof(dfloat));
  parsed->vort_z = (dfloat *) calloc(parsed->points,sizeof(dfloat));

  for(int i = 0; i < parsed->points; i++) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%lf %lf %lf",parsed->vort_x + i,parsed->vort_y+i,parsed->vort_z+i);
  }

  //skip to start of radial vorticity field
  while(!strstr(buf,"RadialVorticity")) fgets(buf,BUFSIZ,fp);  

  parsed->rad_vort = (dfloat *) calloc(parsed->points,sizeof(dfloat));

  for(int i = 0; i < parsed->points; i++) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%lf",parsed->rad_vort+i);
  }

  //skip to start of connectivity field
  while(!strstr(buf,"connectivity")) fgets(buf,BUFSIZ,fp);  

  parsed->connectivity = (iint *) calloc(parsed->elements*4,sizeof(iint));

  for(int i = 0; i < parsed->elements; i++) {
    fgets(buf,BUFSIZ,fp);
    sscanf(buf,"%d %d %d %d",parsed->connectivity+4*i+0,parsed->connectivity+4*i+1,parsed->connectivity+4*i+2,parsed->connectivity+4*i+3);
  }
  
  fclose(fp);

}
