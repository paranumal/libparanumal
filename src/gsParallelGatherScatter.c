/* compile with C compiler (not C++) */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "gslib.h"

void gsParallelGatherScatter(void *gsh, void *v, const char *type, const char *op){

  /* need gs_float or gs_double */
  if(!strcmp(type, "float")){
    //    printf("performing string gs on %s\n", type);
    if(!strcmp(op, "add"))
      gs(v, gs_float, gs_add, 0, gsh, 0);
    else if(!strcmp(op, "min"))
      gs(v, gs_float, gs_min, 0, gsh, 0);
    else if(!strcmp(op, "max"))
      gs(v, gs_float, gs_max, 0, gsh, 0);
  }
  
  if(!strcmp(type, "double")){
    //    printf("performing double gs on %s\n", type);
    if(!strcmp(op, "add"))
      gs(v, gs_double, gs_add, 0, gsh, 0);
    else if(!strcmp(op, "min"))
      gs(v, gs_double, gs_min, 0, gsh, 0);
    else if(!strcmp(op, "max"))
      gs(v, gs_double, gs_max, 0, gsh, 0);
  }

  if(!strcmp(type, "int")){
    //    printf("performing int gs\n");
    if(!strcmp(op, "add"))
      gs(v, gs_int, gs_add, 0, gsh, 0);
    else if(!strcmp(op, "min"))
      gs(v, gs_int, gs_min, 0, gsh, 0);
    else if(!strcmp(op, "max"))
      gs(v, gs_int, gs_max, 0, gsh, 0);
  }
  
}

void gsVecParallelGatherScatter(void *gsh, void *v, const int k, const char *type, const char *op){

  /* need gs_float or gs_double */
  if(!strcmp(type, "float")){
    //    printf("performing string gs on %s\n", type);
    if(!strcmp(op, "add"))
      gs_vec(v, k, gs_float, gs_add, 0, gsh, 0);
    else if(!strcmp(op, "min"))
      gs_vec(v, k, gs_float, gs_min, 0, gsh, 0);
    else if(!strcmp(op, "max"))
      gs_vec(v, k, gs_float, gs_max, 0, gsh, 0);
  }
  
  if(!strcmp(type, "double")){
    //    printf("performing double gs_vec on %s\n", type);
    if(!strcmp(op, "add"))
      gs_vec(v, k, gs_double, gs_add, 0, gsh, 0);
    else if(!strcmp(op, "min"))
      gs_vec(v, k, gs_double, gs_min, 0, gsh, 0);
    else if(!strcmp(op, "max"))
      gs_vec(v, k, gs_double, gs_max, 0, gsh, 0);
  }

  if(!strcmp(type, "int")){
    //    printf("performing int gs_vec\n");
    if(!strcmp(op, "add"))
      gs_vec(v, k, gs_int, gs_add, 0, gsh, 0);
    else if(!strcmp(op, "min"))
      gs_vec(v, k, gs_int, gs_min, 0, gsh, 0);
    else if(!strcmp(op, "max"))
      gs_vec(v, k, gs_int, gs_max, 0, gsh, 0);
  }
  
}
