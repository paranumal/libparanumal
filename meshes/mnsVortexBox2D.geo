res = DefineNumber[0.015];
xmin = DefineNumber[-0.5];
xmax = DefineNumber[ 0.5];
ymin = DefineNumber[-0.5];
ymax = DefineNumber[ 0.5];


Point(1) = {xmin, ymin, 0, res};
Point(2) = {xmin, ymax, 0, res};
Point(3) = {xmax, ymax, 0, res};
Point(4) = {xmax, ymin, 0, res};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(9) = {1, 2, 3, 4};
Plane Surface(9) = {9};
Physical Surface("Domain",9) = {9};
Physical Line("Slip",4) = {1, 2, 3, 4};
