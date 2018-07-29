/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

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
