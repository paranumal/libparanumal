#include "afterparser.h"

int main() {
  data *d = (data *) calloc(1,sizeof(data));

  parse_vtu("./norms_0000_0000.vtu",d);

  send_vtu("./extracted.vtu",d);

  return 0;
}
