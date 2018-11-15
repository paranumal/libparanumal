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


#include "autoTester.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  if(argc!=3){

    printf("usage: ./autoTester executable setupTemplateFile \n");
    exit(-1);
  }
  
  char *executable = strdup(argv[1]);
  setupAide setup(argv[2]);
  char cmd1[BUFSIZ];
  char cmd2[BUFSIZ];
  
  //  sprintf(cmd1, "mpiexec -np 1 %s setup.rc", executable);
  sprintf(cmd1, "%s setup.rc", executable);
  sprintf(cmd2, "mpiexec -np 2 %s setup.rc", executable);
  
  vector <string> &keyword = setup.getKeyword();
  vector <string> &data = setup.getData();

  int Nkeywords = data.size();

  int *counts = (int*) calloc(Nkeywords,sizeof(int));

  int maxNoptions = 100;
  char ***options = (char***) calloc(Nkeywords, sizeof(char**));
  for(int key=0;key<Nkeywords;++key){
    options[key] = (char**) calloc(maxNoptions, sizeof(char*));
  }
  
  for(int key=0;key<Nkeywords;++key){
    printf("%s: %s\n", keyword[key].c_str(), data[key].c_str());

    char *keyStr = strdup(keyword[key].c_str());

    // count number of entries
    char *valStr = strdup(data[key].c_str());
    char *scanStr = strtok(valStr, ",");

    printf("scanStr %s\n", scanStr);
    
    while(scanStr != NULL){

      printf("%s\n", scanStr);

      options[key][counts[key]] = strdup(scanStr);
      
      ++counts[key];
      
      scanStr = strtok(NULL, ",");
    }
  }


  int Ntests = 1;

  for(int key=0;key<Nkeywords;++key){
    Ntests *= counts[key];
  }

  printf("Ntests = %d\n", Ntests);

  int *signature = (int*) calloc(Nkeywords,sizeof(int));

  for(int test=0;test<Ntests;++test){

    int id = 0;

    ++signature[id];
    while(signature[id] == counts[id] && id<Nkeywords){
      signature[id] = 0;
      if(id<Nkeywords-1)
	++signature[id+1];
      ++id;
    }
    fprintf(stderr, "global: ");
    for(id=0;id<Nkeywords;++id){
      fprintf(stderr,"%s,", options[id][signature[id]]);
    }
    fprintf(stderr,"\n");

    FILE *setupFile = fopen("setup.rc", "w");

    for(id=0;id<Nkeywords;++id){
      fprintf(setupFile, "[%s]\n", keyword[id].c_str());
      fprintf(setupFile, "%s\n", options[id][signature[id]]);
    }
    fclose(setupFile);

    fprintf(stderr,"MPI(1): \n");
    system(cmd1);

    //    fprintf(stderr,"MPI(2): \n");
    //    system(cmd2);
    
  }
  
  
  
  MPI_Finalize();
  
  return 1;
}
