#include "afterparser.h"

int main() {
  data *ref = (data *) calloc(1,sizeof(data));
  data *test = (data *) calloc(1,sizeof(data));

  char *test_string = (char *) calloc(300,sizeof(char));
  char *ref_string = (char *) calloc(300,sizeof(char));
  char *out_string = (char *) calloc(300,sizeof(char));
  
  for(int i = 0; i < 63; ++i) {
    
    sprintf(test_string,"/scratch/stimmel/filter/norms_0000_%04d.vtu",i);
    sprintf(ref_string,"/scratch/stimmel/mrab/norms_0000_%04d.vtu",4*(i+1) - 1);
    
    parse_vtu(ref_string,ref);
    parse_vtu(test_string,test);
	    
    dfloat mismatch = compute_l2(ref,test);

    printf("mismatch at test frame %d: %.6lf\n",i,mismatch);
    fflush(stdout);

    sprintf(out_string,"/scratch/stimmel/difference/diff_0000_%04d.vtu",i);
    
    send_vtu(out_string,test);
  }

  return 0;
}
