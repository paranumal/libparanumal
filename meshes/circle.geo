// cl__1 = 0.5;
// Point(1) = {-2,-2, 0, cl__1};
// Point(2) = { 2,-2, 0, cl__1};
// Point(3) = { 2, 2, 0, cl__1};
// Point(4) = {-2, 2, 0, cl__1};

// Line(1) = {1, 2};
// Line(2) = {2, 3};
// Line(3) = {3, 4};
// Line(4) = {4, 1};

// Point(7) = { 0.5, 0.0, 0, cl__1};
// Point(8) = {-0.5, 0.0, 0, cl__1};
// Point(9) = { 0.0, 0.0, 0, cl__1};

// Circle(5) = {8, 9, 7};
// Circle(6) = {7, 9, 8};

// Line Loop(1) = {6, 5};
// Line Loop(2) = {2, 3, 4, 1};

// Plane Surface(1) = {1, 2};
// Physical Surface("Domain", 9) = {1};
// Physical Line("Inflow", 1) = {2, 3, 4, 1};



// cl__1 = 0.3;
cl__1 = 0.25;
Point(1) = {-2,-2, 0, cl__1};
Point(2) = { 2,-2, 0, cl__1};
Point(3) = { 2, 2, 0, cl__1};
Point(4) = {-2, 2, 0, cl__1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {2, 3, 4, 1};

Plane Surface(1) = {1};
Physical Surface("Domain", 9) = {1};
Physical Line("Inflow", 1) = {2, 3, 4, 1};