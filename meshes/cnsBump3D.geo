r1 = DefineNumber[0.10];
r2 = DefineNumber[0.05];
Point(1) = {0, 0, 0, r1};
Point(2) = {1, 0, 0, r2};
Point(3) = {1.5, -1.2, 0, r2};
Point(4) = {2.0, 0.0, 0, r2};
Point(5) = {3.0, 0.0, 0, r1};
Point(6) = {3.0, 1.0, 0, r1};
Point(7) = {0.0, 1.0, 0, r1};
//+
Line(1) = {1, 2};
Circle(2) = {2, 3, 4};
Line(3) = {4, 5};
Line(4) = {5, 6};
Line(5) = {6, 7};
Line(6) = {7, 1};

Line Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
//+
Extrude {0, 0, 1} {Surface{1}; }
Surface Loop(1) = {38, 17, 1, 21, 25, 29, 33, 37};
