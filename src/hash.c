#include <math.h>

unsigned int hash(const unsigned int value) {
  
  const int p[8] = {
    102679, 102701, 102761, 102763,
    102769, 102793, 102797, 102811
  };
  int h[8] = {
    101527, 101531, 101533, 101537,
    101561, 101573, 101581, 101599
  };
    
  const char *ptr = (char*) &value;

  for (int i = 0; i < sizeof(unsigned int); ++i) {
    for (int j = 0; j < 8; ++j) {
      h[j] = (h[j] * p[j])^ptr[i];
    }
  }

  unsigned int ret = h[0];
  for (int i = 1; i < 8; ++i) {
    ret ^= h[i];
  }

  return h[0];
}
