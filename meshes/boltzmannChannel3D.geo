cl__1 = 1;
s = 0.5;
Point(1) = {0,0.5, -0.5, s};
Point(2) = {0,0.5, 0.5, s};
Point(3) = {0, -0.5, 0.5, s};
Point(4) = {0, -0.5, -0.5, s};
Point(5) = {10,0.5, -0.5, s};
Point(6) = {10,0.5, 0.5, s};
Point(7) = {10, -0.5, 0.5, s};
Point(8) = {10, -0.5, -0.5, s};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};

Line(9)  = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Line Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Line Loop(3) = {9, 5, -10, -1};
Plane Surface(3) = {3};

Line Loop(4) = {10, 6, -11, -2};
Plane Surface(4) = {4};

Line Loop(5) = {11, 7, -12, -3};
Plane Surface(5) = {5};

Line Loop(6) = {9, -8, -12, 4};
Plane Surface(6) = {6};

Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

Physical Volume(1) = {1};

Physical Surface("Wall",1) = {3,4,5,6};
Physical Surface("Inflow",2) = {1};
Physical Surface("Outflow",3) = {2};