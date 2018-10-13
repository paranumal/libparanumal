// adapted from https://www.purplealienplanet.com/node/23
#include "types.h"

/* A simple ray tracer */
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h> /* Needed for boolean datatype */
#include <math.h>
#include <omp.h>
#include <cuda.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

/* Width and height of out image */
#define WIDTH  2048
#define HEIGHT 1152
#define DEPTH  2000

#define BOXSIZE 2048

#define TRIANGLE 1
#define SPHERE   2
#define CONE     3
#define DISK     4
#define CYLINDER 5
#define RECTANGLE 6
#define IMAGE     7
#define ELLIPSOID 8

typedef struct {
  // uses bitfields (https://www.geeksforgeeks.org/bit-fields-c/)
  unsigned int refractor: 1; // 1 = launches refraction ray
  unsigned int reflector: 1; // 1 = launches reflection ray
  unsigned int emitter:   1; // 1 = emits light but does not launch forward ray
}info_t;


/* The vector structure */
typedef struct{
  dfloat x,y,z;
}vector_t;

/* The sphere */
typedef struct{
  vector_t pos;
  dfloat  radius;
  vector_t velocity;
  vector_t newVelocity;
  vector_t force;
}sphere_t; 

/* The ellipsoid */
typedef struct{
  vector_t pos;
  vector_t invRadii;
  vector_t velocity;
  vector_t newVelocity;
  vector_t force;
}ellipsoid_t; 


/* The triangle */
typedef struct{
  vector_t vertices[3];
  dfloat   q[3];
}triangle_t;

/* The rectangle */
typedef struct{
  vector_t center;
  vector_t axis[2];
  dfloat   length[2]; // extent of rectangle in axis directions
  int       NI, NJ;
  unsigned char *image;
}rectangle_t;


/* The cone */
typedef struct{
  vector_t vertex; // apex
  vector_t axis;   // axis vector from apex into cone
  dfloat radius;   // base radius
  dfloat height;   // height of apex above base
}cone_t;

/* The disk */
typedef struct{
  vector_t center;
  vector_t normal;
  dfloat   radius;
}disk_t;

/* The cylinder */
typedef struct{
  vector_t center; // center of base disk face
  vector_t axis;   // axis vector from apex into cone
  dfloat radius;   // base radius
  dfloat height;   // height of apex above base
}cylinder_t;


typedef struct{
  dfloat xmin, xmax;
  dfloat ymin, ymax;
  dfloat zmin, zmax;
  int imin, imax;
  int jmin, jmax;
  int kmin, kmax;
}bbox_t;


/* union of shapes */
typedef struct {
  int id;

  int type;
  
  union {
    sphere_t   sphere;
    triangle_t triangle;
    cone_t     cone;
    disk_t     disk;
    cylinder_t cylinder;
    rectangle_t rectangle;
    ellipsoid_t ellipsoid;
  };

  int material;  

  bbox_t bbox;
  
}shape_t;

/* The ray */
typedef struct{
  vector_t start;
  vector_t dir;
  vector_t invDir;
  int level; // which level of the recursion launched this ray
  dfloat   coef;
}ray_t;


/* Colour */
typedef struct{
  dfloat red, green, blue;
}colour_t;

/* Material Definition */
typedef struct{
  colour_t diffuse;
  dfloat reflection;
  dfloat refraction; // transmission coefficient
  dfloat eta;
  info_t info;
}material_t;

/* Lightsource definition */
typedef struct{
  vector_t pos;
  colour_t intensity;
}light_t;

/* sensor configuration */
typedef struct{
  vector_t eyeX;  // model world coordinates of eye
  vector_t Idir; // screen horizontaal direction in model world (unit)
  vector_t Jdir; // screen vertical direction in model world (unit)
  vector_t normal;  // normal through the screen in model world (unit)
  dfloat   Ilength; // length of screen horizontal direction in model world
  dfloat   Jlength; // length of screen vertical direction in model world
  dfloat   offset;
  colour_t bg;      // background color

  dfloat   focalPlaneOffset; // normal distance between sensor and focus plane
  vector_t lensC;       // center of thin lens
  
}sensor_t;


typedef struct{
  int Nshapes;
  shape_t *shapes;
}cell_t;

typedef struct{
  int NI; // number of cells in x direction
  int NJ; // number of cells in y direction
  int NK; // number of cells in z direction
  dfloat xmin;
  dfloat xmax;
  dfloat ymin;
  dfloat ymax;
  dfloat zmin;
  dfloat zmax;
  dfloat dx;
  dfloat dy;
  dfloat dz;
  dfloat invdx;
  dfloat invdy;
  dfloat invdz;

  int      boxCount;
  int     *boxOffsets;
  int *boxContents;
  bbox_t  *bboxes;
  
  int     *c_boxStarts;
  int *c_boxContents;
  bbox_t  *c_bboxes;
}grid_t;

void saveppm(char *filename, unsigned char *img, int width, int height);



__host__ __device__  vector_t vectorCreate(dfloat x, dfloat y, dfloat z);
__host__ __device__  vector_t vectorSub(const vector_t v1, const vector_t v2);
__host__ __device__  dfloat   vectorDot(const vector_t v1, const vector_t v2);
__host__ __device__  vector_t vectorCrossProduct(const vector_t v1, const vector_t v2);
__host__ __device__  vector_t vectorScale(const dfloat c, const vector_t v);
__host__ __device__  vector_t vectorAdd(const vector_t v1, const vector_t v2);
__host__ __device__  dfloat   vectorTripleProduct(const vector_t a, const vector_t b, const vector_t c);
__host__ __device__  vector_t vectorNormalize(const vector_t a);
__host__ __device__  dfloat   vectorNorm(const vector_t a);

typedef struct{

  int Ntriangles;
  
  int Nmaterials;
  material_t *materials;

  int Nshapes;
  shape_t *shapes;
  
  int Nlights;
  light_t *lights;

  grid_t *grid;

  dfloat *randomNumbers;
  dfloat *c_randomNumbers;
  
  unsigned char *c_img;
  shape_t    *c_shapes;
  material_t *c_materials;
  light_t    *c_lights;
  
} scene_t;

scene_t *sceneSetup(int  plotNelements,
		    int  plotField,
		    int *plotEToV,
		    dfloat *plotx,
		    dfloat *ploty,
		    dfloat *plotz,
		    dfloat *plotq);

void render(const scene_t *scene,
	    const dfloat costheta,
	    const dfloat sintheta,
	    unsigned char *img);




void saveImage(char *filename, unsigned char *img, int W, int H);

void simpleRayTracer(int     plotNelements,
		     dfloat *plotx,
		     dfloat *ploty,
		     dfloat *plotz,
		     dfloat *plotq,
		     const char	 *fileBaseName,
		     const int fileIndex);

