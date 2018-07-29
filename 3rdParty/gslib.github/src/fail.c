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

#include <stdio.h>  /* sprintf, vfprintf, stdout */
#include <stdarg.h> /* va_list, va_start, ... */
#include <stdlib.h> /* exit */
#include <string.h> /* memcpy, and str* functions in comm_fail */
#include "name.h"
#include "fail.h"
#include "types.h"
#include "comm.h"

#define nek_exitt FORTRAN_UNPREFIXED(exitt,EXITT)
void die(int status)
{
#ifdef NO_NEK_EXITT
  if(comm_gbl_id==0) exit(status); else for(;;) ;
#else
  //*nek_exitt();
  exit(1);
#endif  
}

void vdiagnostic(const char *prefix, const char *file, unsigned line,
                 const char *fmt, va_list ap)
{
  static char buf[2048]; int n,na,i=0;
  sprintf(buf,"%s(proc %04d, %s:%d): ",prefix,(int)comm_gbl_id,file,line);
  vsprintf(buf+strlen(buf),fmt,ap);
  strcat(buf,"\n");
  n=strlen(buf);
  while(n && (na=fwrite(buf+i,1,n,stdout))) n-=na, i+=na;
  fflush(stdout);
}

void diagnostic(const char *prefix, const char *file, unsigned line,
                const char *fmt, ...)
{
  va_list ap; va_start(ap,fmt);
  vdiagnostic(prefix,file,line,fmt,ap);
  va_end(ap);
}

void vfail(int status, const char *file, unsigned line,
           const char *fmt, va_list ap)
{
  vdiagnostic("ERROR ",file,line,fmt,ap);
  die(status);
}

void fail(int status, const char *file, unsigned line,
          const char *fmt, ...)
{
  va_list ap; va_start(ap,fmt);
  vfail(status,file,line,fmt,ap);
  va_end(ap);
}
