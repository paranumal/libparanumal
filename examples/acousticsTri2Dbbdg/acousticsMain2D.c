#include "acoustics2D.h"

int main(int argc, char **argv){

    // start up MPI
    MPI_Init(&argc, &argv);

    if(argc!=3){
        // to run cavity test case with degree N elements
        printf("usage: ./main meshes/cavityH005.msh N\n");
        exit(-1);
    }

    // int specify polynomial degree 
    int N = atoi(argv[2]);

    // set up mesh stuff
    mesh2D *mesh = meshSetupTri2D(argv[1], N);

    // set up acoustics stuff
    acousticsSetup2D(mesh);

    // run
    //acousticsRun2Dbbdg(mesh);
    acousticsOccaRun2Dbbdg(mesh);

    // close down MPI
    MPI_Finalize();

    exit(0);
    return 0;
}
