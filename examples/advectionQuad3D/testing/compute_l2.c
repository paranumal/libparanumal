#include "afterparser.h"

dfloat compute_l2(data *ref, data *test) {

  //allocate error terms to test code
  test->mismatch = (dfloat *) calloc(ref->points,sizeof(dfloat));
  test->difference = (dfloat *) calloc(ref->points,sizeof(dfloat));
  
  dfloat mismatch = 0;
  for(int p = 0; p < ref->points; ++p) {
    dfloat pw_mismatch = 0;
    dfloat pw_diff = 0;
    for(int fld = 0; fld < 10; ++fld) {
      pw_mismatch += ref->jacobian[p]*(test->q[p*10+fld] - ref->q[p*10+fld])*(test->q[p*10+fld] - ref->q[p*10+fld]);
      pw_diff += (test->q[p*10+fld] - ref->q[p*10+fld])*(test->q[p*10+fld] - ref->q[p*10+fld]);
    }
    mismatch += pw_mismatch;
    test->mismatch[p] = pw_mismatch;
    test->difference[p] = sqrt(pw_mismatch); 
  }
  return sqrt(mismatch);
}
