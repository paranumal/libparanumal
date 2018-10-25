//+
SetFactory("OpenCASCADE");
Sphere(1) = {0, 0, 0, 1, -Pi/2, Pi/2, 2*Pi};
//+
Recombine Surface {1};
