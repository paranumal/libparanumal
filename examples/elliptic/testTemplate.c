
#include "elliptic.h"

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);
  
  setupAide setup(argv[1]);

  matrix <string> &keyword = setup.getKeyword();
  matrix <string> &data = setup.getData();

  int Nkeywords = data.size();

  int *counts = (int*) calloc(Nkeywords,sizeof(int));

  int maxNoptions = 100;
  char ***options = (char***) calloc(Nkeywords, sizeof(char**));
  for(int key=1;key<=Nkeywords;++key){
    options[key-1] = (char**) calloc(maxNoptions, sizeof(char*));
  }
  
  for(int key=1;key<=Nkeywords;++key){
    printf("%s: %s\n", keyword[key].c_str(), data[key].c_str());

    char *keyStr = strdup(keyword[key].c_str());

    // count number of entries
    char *valStr = strdup(data[key].c_str());
    char *scanStr = strtok(valStr, ",");
    
    while(scanStr != NULL){

      printf("%s\n", scanStr);

      options[key-1][counts[key-1]] = strdup(scanStr);
      
      ++counts[key-1];
      
      scanStr = strtok(NULL, ",");
    }
  }


  int Ntests = 1;

  for(int key=1;key<=Nkeywords;++key){
    Ntests *= counts[key-1];
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
      fprintf(setupFile, "[%s]\n", keyword[id+1].c_str());
      fprintf(setupFile, "%s\n", options[id][signature[id]]);
    }
    fclose(setupFile);

    fprintf(stderr,"MPI(1): \n");
    system("mpiexec -n 1 ./ellipticMain setup.rc");

    fprintf(stderr,"MPI(2): \n");
    system("mpiexec -n 2 ./ellipticMain setup.rc");
    
  }
  
  
  
  MPI_Finalize();
  
  return 1;
}
