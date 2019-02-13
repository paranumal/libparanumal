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

constant int p8est_face_edges[6][4] = {{4, 6, 8, 10}, {5, 7, 9, 11},
                                       {0, 2, 8, 9},  {1, 3, 10, 11},
                                       {0, 1, 4, 5},  {2, 3, 6, 7}};
constant int p8est_corner_faces[8][3] = {{0, 2, 4}, {1, 2, 4}, {0, 3, 4},
                                         {1, 3, 4}, {0, 2, 5}, {1, 2, 5},
                                         {0, 3, 5}, {1, 3, 5}};
constant int p8est_corner_edges[8][3] = {{0, 4, 8},  {0, 5, 9}, {1, 4, 10},
                                         {1, 5, 11}, {2, 6, 8}, {2, 7, 9},
                                         {3, 6, 10}, {3, 7, 11}};
constant int p8est_corner_face_corners[8][6] = {
    {0, -1, 0, -1, 0, -1}, {-1, 0, 1, -1, 1, -1}, {1, -1, -1, 0, 2, -1},
    {-1, 1, -1, 1, 3, -1}, {2, -1, 2, -1, -1, 0}, {-1, 2, 3, -1, -1, 1},
    {3, -1, -1, 2, -1, 2}, {-1, 3, -1, 3, -1, 3}};
