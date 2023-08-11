// declare stubbs for missing MPI logging routines.

#pragma once

#ifndef _MSC_VER

#include <string>
#include <cstdio>

void  nnLOG(int n, const char* format_str, ...);
void  nnLOG(const std::string& msg, int n);
char* nnSTR(const char* fmt, ...);

#define nnMSG nnLOG
#define nnTRC nnLOG

// NBN: global g_procid = comm.rank;
extern int g_procid;

#endif
