#include "mesh.h"
#include <stdlib.h>

int isHigher(const void *a, const void *b){

  int *pta = (int*) a;
  int *ptb = (int*) b;

  if(*pta < *ptb) return -1;
  if(*pta > *ptb) return +1;

  return 0;
}

int isLower(const void *a, const void *b){

  int *pta = (int*) a;
  int *ptb = (int*) b;

  if(*pta > *ptb) return -1;
  if(*pta < *ptb) return +1;

  return 0;
}

void mysort(int *data, int N, const char *order){

  if(strstr(order, "ascend")){
    qsort(data, N, sizeof(int), isHigher);
  }
  else{
    qsort(data, N, sizeof(int), isLower);
  }

}
