r = DefineNumber[0.1];
Point(1) = {0, 0, 0, r};
Point(2) = {1, 0, 0, r};
Point(3) = {1, 1, 0, r};
Point(4) = {0, 1, 0, r};
Line(1) = {4, 3};
Line(2) = {2, 3};
Line(3) = {4, 1};
Line(4) = {1, 2};
Line Loop(6) = {3, 4, 2, -1};
Plane Surface(6) = {6};
Recombine Surface {6};
Extrude {0, 0, 1} {
  Surface{6}; Layers{1/r}; Recombine;
}
Physical Surface("Inflow",2) = {28, 6, 27, 19, 15, 23};
Physical Volume("Domain",9) = {1};
