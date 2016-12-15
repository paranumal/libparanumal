#include "mesh3D.h"
#include <stdlib.h>

int isHigher(const void *a, const void *b){

  iint *pta = (iint*) a;
  iint *ptb = (iint*) b;

  if(*pta < *ptb) return -1;
  if(*pta > *ptb) return +1;

  return 0;
}

int isLower(const void *a, const void *b){

  iint *pta = (iint*) a;
  iint *ptb = (iint*) b;

  if(*pta > *ptb) return -1;
  if(*pta < *ptb) return +1;

  return 0;
}

void mysort(iint *data, iint N, const char *order){

  if(strstr(order, "ascend")){
    qsort(data, N, sizeof(iint), isHigher);
  }
  else{
    qsort(data, N, sizeof(iint), isLower);
  }

}
