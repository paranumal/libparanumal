r = DefineNumber[ 0.25];
Point(1) = {0.0, -0.5, -0.5, r};
Point(2) = {5.0, -0.5, -0.5, r};
Point(3) = {5.0,  0.5, -0.5, r};
Point(4) = {0.0,  0.5, -0.5, r};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Line Loop(5) = {1, 2, 3, 4};
//+
Plane Surface(6) = {5};
Recombine Surface {6};
Extrude {0, 0, 1} {
  Surface{6}; Layers{1/r}; Recombine;
}
Physical Surface("Wall",1) = {28, 6, 23, 15};
Physical Surface("Inflow",2) = {27};
Physical Surface("Outflow",3) = {19};
Physical Volume("Domain",9) = {1};
