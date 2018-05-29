cl__1 = 1;

r_0 = DefineNumber[1.0];
r_1 = DefineNumber[1.0];

Nc  = DefineNumber[30];
Np  = DefineNumber[5];
Nb  = DefineNumber[10];

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

Point(21) = {-xmin, ymax, 0, r_1};
Point(22) = {-xmin, ymin, 0, r_1};

Point(23) = {-xmin, ypmlmax, 0, r_1};
Point(24) = {-xmin, ypmlmin, 0, r_1};

Line(1) = {17, 18};
Line(2) = {18, 19};
Line(3) = {19, 20};
Line(4) = {20, 17};
Line(5) = {19, 21};
Line(6) = {21, 4};
Line(7) = {4, 20};
Line(8) = {17, 1};
Line(9) = {1, 4};
Line(10) = {1, 22};
Line(11) = {22, 18};
Line(12) = {22, 21};
Line(13) = {21, 3};
Line(14) = {3, 2};
Line(15) = {2, 22};
Line(16) = {3, 14};
Line(17) = {14, 12};
Line(18) = {12, 2};
Line(19) = {12, 6};
Line(20) = {6, 11};
Line(21) = {11, 2};
Line(22) = {11, 24};
Line(23) = {24, 22};
Line(24) = {24, 9};
Line(25) = {9, 1};
Line(26) = {9, 5};
Line(27) = {5, 10};
Line(28) = {10, 1};
Line(29) = {10, 16};
Line(30) = {4, 16};
Line(31) = {16, 8};
Line(32) = {8, 15};
Line(33) = {15, 4};
Line(34) = {15, 23};
Line(35) = {23, 21};
Line(36) = {23, 13};
Line(37) = {13, 3};
Line(38) = {13, 7};
Line(39) = {7, 14};


Transfinite Line {7, 8, 5, 11} = Nc Using Progression 1;
Transfinite Line {29, 9, 4, 2, 12, 14, 17} = Nc Using Progression 1;
Transfinite Line {26, 28, 30, 32, 38, 16, 18, 20} = Np Using Progression 1;
Transfinite Line {22, 15, 13, 36} = Nb Using Progression 1;
Transfinite Line {24, 10, 1, 3, 6, 34} = Nc Using Progression 1;
Transfinite Line {31, 33, 35, 37, 39} = Np Using Progression 1;
Transfinite Line {27, 25, 23, 21, 19} = Np Using Progression 1;

Line Loop(40) = {3, -7, -6, -5};
Plane Surface(41) = {40};
Line Loop(42) = {2, 5, -12, 11};
Plane Surface(43) = {42};
Line Loop(44) = {11, -1, 8, 10};
Plane Surface(45) = {44};
Line Loop(46) = {8, 9, 7, 4};
Plane Surface(47) = {46};
Line Loop(48) = {28, -25, 26, 27};
Plane Surface(49) = {48};
Line Loop(50) = {24, 25, 10, -23};
Plane Surface(51) = {50};
Line Loop(52) = {23, -15, -21, 22};
Plane Surface(53) = {52};
Line Loop(54) = {18, -21, -20, -19};
Plane Surface(55) = {54};
Line Loop(56) = {14, 15, 12, 13};
Plane Surface(57) = {56};
Line Loop(58) = {14, -18, -17, -16};
Plane Surface(59) = {58};
Line Loop(60) = {39, -16, -37, 38};
Plane Surface(61) = {60};
Line Loop(62) = {13, -37, -36, 35};
Plane Surface(63) = {62};
Line Loop(64) = {6, -33, 34, 35};
Plane Surface(65) = {64};
Line Loop(66) = {33, 30, 31, 32};
Plane Surface(67) = {66};
Line Loop(68) = {9, 30, -29, 28};
Plane Surface(69) = {68};
Transfinite Surface {69} = {1, 10, 16, 4};
Transfinite Surface {67} = {4, 16, 8, 15};
Transfinite Surface {65} = {4, 15, 23, 21};
Transfinite Surface {63} = {21, 23, 13, 3};
Transfinite Surface {61} = {13, 7, 14, 3};
Transfinite Surface {59} = {3, 14, 12, 2};
Transfinite Surface {55} = {12, 6, 11, 2};
Transfinite Surface {53} = {11, 24, 22, 2};
Transfinite Surface {51} = {22, 24, 9, 1};
Transfinite Surface {49} = {1, 9, 5, 10};
Transfinite Surface {57} = {21, 3, 22, 2};
Transfinite Surface {43} = {19, 21, 18, 22};
Transfinite Surface {41} = {20, 4, 19, 21};
Transfinite Surface {47} = {20, 4, 17, 1};
Transfinite Surface {45} = {17, 1, 18, 22};

Recombine Surface {67, 65, 63, 61, 59, 55, 53, 57, 51, 49, 69, 47, 41, 43, 45};

Coherence;


Physical Line("Wall",1) = {3, 2, 1, 4};
Physical Line("Outflow",3) = {29, 31, 32, 34, 36, 38, 39, 17, 20, 19, 22, 24, 26, 27};

Physical Surface("Interior",9) = {47, 41, 43, 45, 57};
Physical Surface("XPML",100) = {69, 59};
Physical Surface("YPML",200) = {51, 53, 63, 65};
Physical Surface("XYPML",300) = {49, 67, 61, 55};
