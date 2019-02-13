#define mfloat dfloat

#define APPEND(x, y) x##y
#define APPEND_EXPAND(x, y) BFAM_APPEND(x, y)

#ifdef p_DFLOAT_FLOAT
#define DFLOAT(x) APPEND(x, f)
#else
#define DFLOAT(x) x
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define ABS(a) (((a) > 0) ? (a) : -(a))

#define ZERO DFLOAT(0.0)
#define ONE DFLOAT(1.0)
#define TWO DFLOAT(2.0)
#define HALF DFLOAT(0.5)

