r = DefineNumber[ 0.15];
Point(1) = {-0.5, -0.5, 0.0, r};
Point(2) = {-0.5, 1.5,  0.0, r};
Point(3) = {1.0, -0.5,  0.0, r};
Point(4) = {1.0, 1.5, 0.0, r};

Line(1) = {2, 4};
Line(2) = {4, 3};
Line(3) = {3, 1};
Line(4) = {1, 2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Extrude {0, 0, 2} { Surface{6};}

Physical Surface("Inflow", 2) = {15, 6, 23, 27, 28};
Physical Surface("Outflow", 3) = {19};
Physical Volume("Domain", 9) = {1};

//+Mesh.RemeshAlgorithm = 1; // automatic
//+Mesh.RemeshParametrization = 7; // conformal finite element
//+Mesh.Algorithm = 6; // Frontal

