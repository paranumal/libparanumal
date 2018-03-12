cl__1 = 0.1;
Point(1) = {-1, -1, -1, cl__1};
Point(2) = {1, -1, -1, cl__1};
Point(3) = {1, 1, -1, cl__1};
Point(4) = {-1, 1, -1, cl__1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(6) = {4, 1, 2, 3};
Plane Surface(6) = {6};
Recombine Surface {6};
Extrude {0, 0, 2} {
  Surface{6}; Layers{2/cl__1}; Recombine;
}
Physical Surface("Inflow", 1) = {28, 23, 6, 15, 19, 27};
Physical Volume("Domain", 9) = {1};
