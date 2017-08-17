cl__1 = 1;

Point(1) = {-1, 0, 0, cl__1};
Point(2) = { 1, 0, 0, cl__1};
Point(3) = { 0, 1.732050808, 0, cl__1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 1};
Line Loop(4) = {1, 2, 3};
Plane Surface(5) = {4};
Physical Surface("Domain", 9) = {5};
Physical Line("Inflow", 1) = {1, 2, 3};
