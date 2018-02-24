cl__1 = 1;

r_0 = DefineNumber[0.25];
r_1 = DefineNumber[0.025];
r_2 = DefineNumber[4.0];

xmin = DefineNumber[-5.4];
xmax = DefineNumber[ 9.4];
ymin = DefineNumber[-5.4];
ymax = DefineNumber[ 5.4];

xpmlmin = DefineNumber[xmin - 2.6];
xpmlmax = DefineNumber[xmax + 2.6]; 
ypmlmin = DefineNumber[ymin - 2.6];
ypmlmax = DefineNumber[ymax + 2.6];

xdmin = DefineNumber[-399.0];
xdmax = DefineNumber[ 399.0]; 
ydmin = DefineNumber[-399.0];
ydmax = DefineNumber[ 399.0];


xcmin = DefineNumber[-0.5];
xcmax = DefineNumber[ 0.5];
ycmin = DefineNumber[-0.5];
ycmax = DefineNumber[ 0.5];

Point(1) = {xmin, ymin, 0, r_0};
Point(2) = {xmax, ymin, 0, r_0};
Point(3) = {xmax, ymax, 0, r_0};
Point(4) = {xmin, ymax, 0, r_0};

Point(5) = {xpmlmin, ypmlmin, 0, r_0};
Point(6) = {xpmlmax, ypmlmin, 0, r_0};
Point(7) = {xpmlmax, ypmlmax, 0, r_0};
Point(8) = {xpmlmin, ypmlmax, 0, r_0};

Point(9) = {xmin, ypmlmin, 0, r_0};
Point(10) = {xpmlmin, ymin, 0, r_0};

Point(11) = {xmax, ypmlmin, 0, r_0};
Point(12) = {xpmlmax, ymin, 0, r_0};

Point(13) = {xmax, ypmlmax, 0, r_0};
Point(14) = {xpmlmax, ymax, 0, r_0};

Point(15) = {xmin, ypmlmax, 0, r_0};
Point(16) = {xpmlmin, ymax, 0, r_0};

Point(17) = {xcmin, ycmin, 0, r_1};
Point(18) = {xcmax, ycmin, 0, r_1};
Point(19) = {xcmax, ycmax, 0, r_1};
Point(20) = {xcmin, ycmax, 0, r_1};

Point(21) = {xdmin, ydmin, 0, r_2};
Point(22) = {xdmax, ydmin, 0, r_2};
Point(23) = {xdmax, ydmax, 0, r_2};
Point(24) = {xdmin, ydmax, 0, r_2};


Line(1) = {17, 18};
Line(2) = {18, 19};
Line(3) = {19, 20};
Line(4) = {20, 17};
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};
Line(9) = {5, 9};
Line(10) = {9, 1};
Line(11) = {1, 10};
Line(12) = {10, 5};
Line(13) = {9, 11};
Line(14) = {2, 11};
Line(15) = {11, 6};
Line(16) = {6, 12};
Line(17) = {12, 2};
Line(18) = {12, 14};
Line(19) = {14, 3};
Line(20) = {14, 7};
Line(21) = {7, 13};
Line(22) = {13, 3};
Line(23) = {13, 15};
Line(24) = {15, 4};
Line(25) = {15, 8};
Line(26) = {8, 16};
Line(27) = {16, 4};
Line(28) = {16, 10};

Line(29) = {21, 24};
Line(30) = {24, 23};
Line(31) = {23, 22};
Line(32) = {22, 21};
Line Loop(33) = {7, 8, 5, 6};
Line Loop(34) = {3, 4, 1, 2};
Plane Surface(35) = {33, 34};
Line Loop(36) = {5, 14, -13, 10};
Plane Surface(37) = {36};
Line Loop(38) = {12, 9, 10, 11};
Plane Surface(39) = {38};
Line Loop(40) = {28, -11, -8, -27};
Plane Surface(41) = {40};
Line Loop(42) = {26, 27, -24, 25};
Plane Surface(43) = {42};
Line Loop(44) = {23, 24, -7, -22};
Plane Surface(45) = {44};
Line Loop(46) = {21, 22, -19, 20};
Plane Surface(47) = {46};
Line Loop(48) = {18, 19, -6, -17};
Plane Surface(49) = {48};
Line Loop(50) = {15, 16, 17, 14};
Plane Surface(51) = {50};
Line Loop(52) = {25, 26, 28, 12, 9, 13, 15, 16, 18, 20, 21, 23};
Line Loop(53) = {32, 29, 30, 31};
Plane Surface(54) = {52, 53};

Physical Surface("Interior",9) = {54, 39, 37, 51, 49, 47, 45, 41, 43, 35};

Physical Line("Inflow",2) = {29, 30, 32};
Physical Line("Outflow",3) = {31};
Physical Line("Wall",1) = {1, 2, 3, 4};

Coherence;
