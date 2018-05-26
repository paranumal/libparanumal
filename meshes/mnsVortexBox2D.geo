res = DefineNumber[0.025];
Point(1) = {0, 0, 0, res};
Point(2) = {0, 1.0, 0, res};
Point(3) = {1.0, 1.0, 0, res};
Point(4) = {1.0, 0, 0, res};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(9) = {1, 2, 3, 4};
Plane Surface(9) = {9};
Physical Surface("Domain",9) = {9};
Physical Line("Slip",4) = {1, 2, 3, 4};
