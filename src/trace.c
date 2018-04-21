#include<trace.hpp>

using std::string;
// GNU specialization. Not portable.

ostream &trace(ostream &stream, int stack_size){

  void *array[stack_size];

  int length = backtrace (array, stack_size);
  char **trace = backtrace_symbols (array, length);

  //  for( size_t i=1; i<length; ++i ){
  //    cout << trace[i] << endl;
  //  }
  //  return stream;

  string binary;

  // start from 1, since the first entry is the traced_exception ctor
  for( size_t i=1; i<length; ++i ){
    char *line = trace[i];
    char *lb = strchr(line, '(');
    char *plus = 0;
    char *rb = 0;
    char *demangled = 0;
    if( lb != 0 ){
      if( binary != string( line, lb-line ) ){
	binary = string( line, lb-line );
	stream << "in " << binary << ":\n";
      }
      *lb = '\0';
      plus = strchr(lb+1,'+');
      if( plus != 0 ){
	*plus = '\0';
	int status;
	demangled = abi::__cxa_demangle(lb+1, 0, 0, &status);
	rb = strchr(plus+1, ')');
      }else{
	rb = strchr(lb+1, ')');
      }
      if( rb != 0 ){
	*rb = '\0';
      }
    }
    if( lb == 0 ){
      stream << line;
    }else{
      if( plus != 0 ){
	stream << "+" << (plus+1);
      }
      if( demangled != 0 ){
	stream << '\t' << demangled;
	free(demangled);
      }else if( lb != 0 && rb > lb ){
	stream << '\t' << (lb+1) << "()";
      }
    }
    stream << '\n';
  }

  return stream;
}
