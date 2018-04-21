#ifndef __TRACE
#define __TRACE

// using backtrace
// http://www.gnu.org/software/libc/manual/html_node/Backtraces.html
// demangling output:
// http://gcc.gnu.org/onlinedocs/libstdc++/manual/ext_demangling.html
// spec: http://www.ib.cnea.gov.ar/~oop/biblio/libstdc++/namespaceabi.html


#include <iostream>
#include <execinfo.h>
#include <stdlib.h>
#include <cstring>
#include <cxxabi.h>

using std::ostream;
//#include <headers2d.hpp>
// GNU specialization. Not portable.
ostream &trace(ostream &stream, int stack_size);

#endif
