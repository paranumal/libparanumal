// define stubbs for missing MPI logging routines.

#ifndef _MSC_VER

#include "stubb_log.hpp"

int g_procid = 0;     // g_procid = comm.rank;
int g_LOG_FLAG = 5;   // manage output

// using standard printf syntax:
void nnLOG(int n, const char* format_str, ...) {

  static char bufL[1024];

  // optionally, only rank 0 prints messages
  bool b = (0 == g_procid);

  if (b && (n <= g_LOG_FLAG)) {
    va_list arglist;
    va_start(arglist, format_str);
    int nUsed = -1;
    nUsed = vsnprintf  (bufL, 1023,       format_str, arglist);
  //nUsed = vsnprintf_s(bufL, 1023, 1000, format_str, arglist);
    assert(nUsed>=0);
    va_end(arglist);

    fprintf(stdout, "%s", bufL);
    fflush (stdout);
  }
}

// handle std::strings
void nnLOG(const std::string& msg, int n) {
  nnLOG(n, "%s", msg.c_str());
}

// Internal formatting of small strings using internal buffer 
char* nnSTR(const char* fmt, ...) {
  va_list ap;
  va_start(ap, fmt);
  static int ncalls = 0;
  const  int repeat_cycle = 8;
  static char s [repeat_cycle][999];
  if (repeat_cycle == ++ncalls) {
    ncalls = 0;
  }
  vsprintf  (s[ncalls],      fmt, ap);
//vsprintf_s(s[ncalls], 900, fmt, ap);
  va_end(ap);

  return s[ncalls];
}

#endif
