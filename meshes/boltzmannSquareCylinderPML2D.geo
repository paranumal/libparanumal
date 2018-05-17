cl__1 = 1;

r_0 = DefineNumber[1.0];
r_1 = DefineNumber[0.25];


xmin = DefineNumber[-5.4];
xmax = DefineNumber[ 9.4];
ymin = DefineNumber[-5.4];
ymax = DefineNumber[ 5.4];

xpmlmin = DefineNumber[xmin - 1*1.6];
xpmlmax = DefineNumber[xmax + 1*1.6]; 
ypmlmin = DefineNumber[ymin - 1*1.6];
ypmlmax = DefineNumber[ymax + 1*1.6];


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





Line Loop(29) = {5, 6, 7, 8};
Line Loop(30) = {1, 2, 3, 4};
Plane Surface(31) = {29, 30};
Line Loop(32) = {9, 10, 11, 12};
Plane Surface(33) = {32};
Line Loop(34) = {10, 5, 14, -13};
Plane Surface(35) = {34};
Line Loop(36) = {15, 16, 17, 14};
Plane Surface(37) = {36};
Line Loop(38) = {18, 19, -6, -17};
Plane Surface(39) = {38};
Line Loop(40) = {20, 21, 22, -19};
Plane Surface(41) = {40};
Line Loop(42) = {22, 7, -24, -23};
Plane Surface(43) = {42};
Line Loop(44) = {25, 26, 27, -24};
Plane Surface(45) = {44};
Line Loop(46) = {28, -11, -8, -27};
Plane Surface(47) = {46};




Physical Line("Wall",1) = {1, 2, 3, 4};
Physical Line("Inflow",2) = {12, 28, 26, 25, 23, 21, 15, 13, 9};
Physical Line("Outflow",3) = {16, 18, 20};



Physical Surface("Interior",9) = {31};
Physical Surface("XPML",100) = {39,47};
Physical Surface("YPML",200) = {35,43};
Physical Surface("XYPML",300) = {33, 37, 41, 45};

Coherence;
