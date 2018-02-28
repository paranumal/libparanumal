#include "mesh.h"
#include <stdlib.h>

int isHigher(const void *a, const void *b){

  hlong *pta = (hlong*) a;
  hlong *ptb = (hlong*) b;

  if(*pta < *ptb) return -1;
  if(*pta > *ptb) return +1;

  return 0;
}

int isLower(const void *a, const void *b){

  hlong *pta = (hlong*) a;
  hlong *ptb = (hlong*) b;

  if(*pta > *ptb) return -1;
  if(*pta < *ptb) return +1;

  return 0;
}

void mysort(hlong *data, int N, const char *order){

  if(strstr(order, "ascend")){
    qsort(data, N, sizeof(hlong), isHigher);
  }
  else{
    qsort(data, N, sizeof(hlong), isLower);
  }

}
