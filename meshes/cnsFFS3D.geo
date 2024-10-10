cl__1 = 0.025;
cl__2 = 0.0075;
Point(1) = {0, 0, 0, cl__1};
Point(2) = {0.6, 0, 0, cl__2};
Point(3) = {0.6, 0.2, 0, cl__2};
Point(4) = {3, 0.2, 0, cl__1};
Point(5) = {3, 1, 0, cl__1};
Point(6) = {0, 1, 0, cl__1};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};
