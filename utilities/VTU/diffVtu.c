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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){

  if(argc!=4){
    printf("usage: ./diffVtu foo1.vtu foo2.vtu foo2m1.vtu\n");
    exit(-1);
  }
    

  FILE *in1 = fopen(argv[1], "r");
  FILE *in2 = fopen(argv[2], "r");
  FILE *out = fopen(argv[3], "w");

  // scan header
  char buf1[BUFSIZ];
  char buf2[BUFSIZ];

  // keep printing lines from file 1 until hit Pressure
  do{
    
    fgets(buf1, BUFSIZ, in1);
    fgets(buf2, BUFSIZ, in2);
    fprintf(out, "%s", buf1);

  }while(!strstr(buf1, "ressure"));

  do{
    float p1,p2;
    fgets(buf1, BUFSIZ, in1);
    fgets(buf2, BUFSIZ, in2);
    if(!strstr(buf1,">")){
      sscanf(buf1, "%f", &p1);
      sscanf(buf2, "%f", &p2);
      fprintf(out, "%f\n", p2-p1);
    }
    else
      fprintf(out, "%s", buf1);
    
  } while(!strstr(buf1, ">"));


  do{
    
    fgets(buf1, BUFSIZ, in1);
    fprintf(out, "%s", buf1);
    
  }while(!strstr(buf1, "VTKFile"));

  fprintf(out, "\n");
  fclose(in1);
  fclose(in2);
  fclose(out);

  return 0;
}
