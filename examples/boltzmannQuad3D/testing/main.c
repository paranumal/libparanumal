#include "afterparser.h"

int main() {
  data *ref = (data *) calloc(1,sizeof(data));
  data *test = (data *) calloc(1,sizeof(data));

  char *test_string = (char *) calloc(300,sizeof(char));
  char *ref_string = (char *) calloc(300,sizeof(char));
  char *out_string = (char *) calloc(300,sizeof(char));
  
  for(int i = 0; i < 18; ++i) {
    
    sprintf(test_string,"/scratch/stimmel/short_filter/norms_0000_%04d.vtu",i);
    sprintf(ref_string,"/scratch/stimmel/short_mrab/norms_0000_%04d.vtu",(i+1) - 1);
    
    parse_vtu(ref_string,ref);
    parse_vtu(test_string,test);
	    
    dfloat mismatch = compute_l2(ref,test);

    printf("mismatch at test frame %d: %.10lf\n",i,mismatch);
    fflush(stdout);

    sprintf(out_string,"/scratch/stimmel/short_diff/diff_0000_%04d.vtu",i);
    
    send_vtu(out_string,test);
    //free most egregious memory offenders
    free(ref->x);
    free(ref->y);
    free(ref->z);
    free(ref->mrab_levels);
    free(ref->q);
    free(ref->density);
    free(ref->jacobian);
    free(ref->vel_x);
    free(ref->vel_y);
    free(ref->vel_z);
    free(ref->vort_x);
    free(ref->vort_y);
    free(ref->vort_z);
    free(ref->rad_vort);
    free(ref->connectivity);

    free(test->mismatch);
    free(test->difference);
    free(test->x);
    free(test->y);
    free(test->z);
    free(test->mrab_levels);
    free(test->q);
    free(test->density);
    free(test->jacobian);
    free(test->vel_x);
    free(test->vel_y);
    free(test->vel_z);
    free(test->vort_x);
    free(test->vort_y);
    free(test->vort_z);
    free(test->rad_vort);
    free(test->connectivity);    
  }

  return 0;
}
