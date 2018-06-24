r0 = 0.010;
r1 = 0.025;
r2 = 0.075;
r3 = 0.15;

ld = 0.5;

// x-bounds
xdmin  = -8.0;
xdmax  = 25.0;
pmlx   =  4.0;
// y-bounds
ydmin  = -8.0;
ydmax  =  8.0;
pmly   =  4.0;

// z-bounds
zall = 0.0;
//------ Domain-------------------------
Point(1) = { xdmin,  ydmin, zall, ld};
Point(2) = { xdmax,  ydmin, zall, r3};
Point(3) = { xdmin,  ydmax, zall, ld};
Point(4) = { xdmax,  ydmax, zall, r3};


th = 1.0;            // height of fence
tw = 0.1;            // width of fence
tt = 5.0;            //legth of tail
//---------------------------------------
Point(5)  = { -tw/2, -th   , zall, r0};  // base front
Point(6)  = {  tw/2, -th   , zall, r0};  // base front
Point(7)  = {  tw/2, -tw/2,  zall, r0};  // base front
Point(8)  = {  tt  , -tw/2,  zall, r0};  // base front
Point(9)  = {  tt  ,  tw/2,  zall, r0};  // base front
Point(10) = {  tw/2,  tw/2,  zall, r0};  // base front
Point(11) = {  tw/2,  th   , zall, r0};  // base front
Point(12) = { -tw/2,  th   , zall, r0};  // base front

th2 = 1.75;           // height of fence
tw2 = 1.0;            // width of fence
tt2 = 8.0;            //legth of tail
//---------------------------------------
Point(13)  = { -tw2/2, -th2   , zall, r1};  // base front
Point(14)  = {  tw2/2, -th2   , zall, r1};  // base front
Point(15)  = {  tw2/2, -tw2/2,  zall, r1};  // base front
Point(16)  = {  tt2  , -tw2/2,  zall, r1};  // base front
Point(17)  = {  tt2  ,  tw2/2,  zall, r1};  // base front
Point(18)  = {  tw2/2,  tw2/2,  zall, r1};  // base front
Point(19)  = {  tw2/2,  th2   , zall, r1};  // base front
Point(20)  = { -tw2/2,  th2   , zall, r1};  // base front


// Create A second smmoth transition region
xr1  = 10.0;
yr1  = 7.0; 

// Point(21)  = { xr1,  yr1   , zall, r2};  // base front
// Point(22)  = { xr1, -yr1   , zall, r2};  // base front

// Dummy points for mesh

Point(41)  = { xr1,   yr1   , zall, r2};  // base front
Point(42)  = { xr1,  -yr1   , zall, r2};  // base front
Point(43)  = { xr1,  ydmax+pmly   , zall, r3};  // base front
Point(44)  = { xr1,  ydmin-pmly   , zall, r3};  // base front
Point(45)  = { xr1,   ydmax   , zall, r3};  // base front
Point(46)  = { xr1,   ydmin   , zall, r3};  // base front

//------ Domain-------------------------
Point(23) = { xdmin-pmlx,  ydmin, zall, ld};
Point(24) = { xdmax+pmlx,  ydmin, zall, r3};
Point(25) = { xdmin-pmlx,  ydmax, zall, ld};
Point(26) = { xdmax+pmlx,  ydmax, zall, r3};
//
Point(27) = { xdmin,  ydmin-pmly, zall, ld};
Point(28) = { xdmax,  ydmin-pmly, zall, r3};
Point(29) = { xdmin,  ydmax+pmly, zall, ld};
Point(30) = { xdmax,  ydmax+pmly, zall, r3};
//
Point(31) = { xdmin-pmlx,  ydmin-pmly, zall, ld};
Point(32) = { xdmax+pmlx,  ydmin-pmly, zall, r3};
Point(33) = { xdmin-pmlx,  ydmax+pmly, zall, ld};
Point(34) = { xdmax+pmlx,  ydmax+pmly, zall, r3};

// // 
// yr2  = 7.0; 
// Point(35)  = { xdmax, yr2, zall, r3};  // base front
// Point(36)  = { xdmax,-yr2, zall, r3};  // base front
// Point(37)  = { xdmax+pmlx, yr2, zall, r3};  // base front
// Point(38)  = { xdmax+pmlx,-yr2, zall, r3};  // base front


yr3 = 2.75;
Point(39)  = { -tw2/2,  yr3   , zall, r1};  // base front
Point(40)  = { -tw2/2, -yr3   , zall, r1};  // base front


// Line(1) = {5, 6};
// Line(2) = {6, 7};
// Line(3) = {7, 8};
// Line(4) = {8, 9};
// Line(5) = {9, 10};
// Line(6) = {10, 11};
// Line(7) = {11, 12};
// Line(8) = {12, 5};
// Line(9) = {13, 14};
// Line(10) = {14, 15};
// Line(11) = {15, 16};
// Line(12) = {16, 17};
// Line(13) = {17, 18};
// Line(14) = {18, 19};
// Line(15) = {19, 20};
// Line(16) = {20, 13};
// Line(17) = {20, 39};
// Line(18) = {39, 21};
// Line(19) = {21, 35};
// Line(20) = {13, 40};
// Line(21) = {40, 22};
// Line(22) = {22, 36};
// Line(23) = {35, 36};
// Line(24) = {35, 4};
// Line(25) = {4, 30};
// Line(26) = {30, 34};
// Line(27) = {34, 26};
// Line(28) = {26, 4};
// Line(29) = {26, 37};
// Line(30) = {37, 35};
// Line(31) = {37, 38};
// Line(32) = {38, 36};
// Line(33) = {36, 2};
// Line(34) = {2, 24};
// Line(35) = {38, 24};
// Line(36) = {24, 32};
// Line(37) = {32, 28};
// Line(38) = {28, 2};
// Line(39) = {4, 3};
// Line(40) = {3, 29};
// Line(41) = {29, 30};
// Line(42) = {3, 25};
// Line(43) = {25, 33};
// Line(44) = {29, 33};
// Line(45) = {25, 23};
// Line(46) = {23, 1};
// Line(47) = {1, 3};
// Line(48) = {23, 31};
// Line(49) = {31, 27};
// Line(50) = {27, 1};
// Line(51) = {1, 2};
// Line(52) = {28, 27};
// Line Loop(53) = {5, 6, 7, 8, 1, 2, 3, 4};
// Line Loop(54) = {14, 15, 16, 9, 10, 11, 12, 13};
// Plane Surface(55) = {53, 54};
// Line Loop(56) = {15, 17, 18, 19, 23, -22, -21, -20, 9, 10, 11, 12, 13, 14};
// Plane Surface(57) = {56};
// Line Loop(58) = {23, -32, -31, 30};
// Plane Surface(59) = {58};
// Line Loop(60) = {24, -28, 29, 30};
// Plane Surface(61) = {60};
// Line Loop(62) = {32, 33, 34, -35};
// Plane Surface(63) = {62};
// Line Loop(64) = {36, 37, 38, 34};
// Plane Surface(65) = {64};
// Line Loop(66) = {25, 26, 27, 28};
// Line Loop(67) = {52, 50, 51, -38};
// Plane Surface(68) = {66};
// Line Loop(69) = {39, 40, 41, -25};
// Plane Surface(70) = {69};
// Line Loop(71) = {44, -43, -42, 40};
// Plane Surface(72) = {71};
// Line Loop(73) = {42, 45, 46, 47};
// Plane Surface(74) = {73};
// Line Loop(75) = {48, 49, 50, -46};
// Plane Surface(76) = {75};
// Plane Surface(77) = {67};
// Line Loop(78) = {24, 39, -47, 51, -33, -22, -21, -20, -16, 17, 18, 19};
// Plane Surface(79) = {78};


//Physical Line("Wall", 1) = {5, 4, 3, 2, 1, 8, 7, 6};
//Physical Line("Slip", 5) = {49, 52, 44, 41, 26,37};
//Physical Line("Outflow", 3) = {27, 29, 31, 35, 36};
//Physical Line("Inflow", 2) = {43, 45, 48};
//Physical Surface("XPML",100) = {74, 61, 59, 63};
//Physical Surface("YPML",200) = {70, 77};
//Physical Surface("XYPML",300) = {72, 68, 65, 76};
//Physical Surface("Interior",9) = {79, 57, 55};


Line(1) = {5, 6};
Line(2) = {6, 7};
Line(3) = {7, 8};
Line(4) = {8, 9};
Line(5) = {9, 10};
Line(6) = {10, 11};
Line(7) = {11, 12};
Line(8) = {12, 5};
Line(9) = {13, 14};
Line(10) = {14, 15};
Line(11) = {15, 16};
Line(12) = {16, 17};
Line(13) = {17, 18};
Line(14) = {18, 19};
Line(15) = {19, 20};
Line(16) = {20, 13};
Line(17) = {13, 40};
Line(18) = {40, 42};
Line(19) = {42, 46};
Line(20) = {20, 39};
Line(21) = {39, 41};
Line(22) = {41, 45};
Line(23) = {45, 4};
Line(24) = {4, 2};
Line(25) = {2, 46};
Line(26) = {46, 44};
Line(27) = {44, 28};
Line(28) = {28, 32};
Line(29) = {32, 24};
Line(30) = {24, 26};
Line(31) = {26, 34};
Line(32) = {34, 30};
Line(33) = {30, 4};
Line(34) = {4, 26};
Line(35) = {24, 2};
Line(36) = {2, 28};
Line(37) = {43, 30};
Line(38) = {43, 45};
Line(39) = {45, 3};
Line(40) = {3, 29};
Line(41) = {29, 43};
Line(43) = {29, 33};
Line(44) = {33, 25};
Line(45) = {25, 3};
Line(46) = {3, 1};
Line(47) = {1, 23};
Line(48) = {23, 25};
Line(50) = {23, 31};
Line(51) = {31, 27};
Line(52) = {27, 1};
Line(53) = {27, 44};
Line(54) = {46, 1};
Line Loop(55) = {5, 6, 7, 8, 1, 2, 3, 4};
Line Loop(56) = {15, 16, 9, 10, 11, 12, 13, 14};
Plane Surface(57) = {55, 56};
Line Loop(58) = {18, 19, -25, -24, -23, -22, -21, -20, 16, 17};
Line Loop(59) = {23, 24, 25, -19, -18, -17, 9, 10, 11, 12, 13, 14, 15, 20, 21, 22};
Plane Surface(60) = {59};
Line Loop(61) = {39, 46, -54, -19, -18, -17, -16, 20, 21, 22};
Plane Surface(62) = {61};
Line Loop(63) = {37, 33, -23, -38};
Plane Surface(64) = {63};
Line Loop(65) = {32, 33, 34, 31};
Plane Surface(66) = {65};
Line Loop(67) = {30, -34, 24, -35};
Plane Surface(68) = {67};
Line Loop(69) = {29, 35, 36, 28};
Plane Surface(70) = {69};
Line Loop(71) = {27, -36, 25, 26};
Plane Surface(72) = {71};
Line Loop(73) = {54, -52, 53, -26};
Plane Surface(74) = {73};
Line Loop(75) = {47, 50, 51, 52};
Plane Surface(76) = {75};
Line Loop(77) = {46, 47, 48, 45};
Plane Surface(78) = {77};
Line Loop(79) = {40, 43, 44, 45};
Plane Surface(80) = {79};
Line Loop(81) = {41, 38, 39, 40};
Plane Surface(82) = {81};

Physical Line("Wall", 1) = {5, 4, 3, 2, 1, 8, 7, 6};
Physical Line("Slip", 5) = {43, 41, 37, 32, 28, 27, 53, 51};
Physical Line("Outflow", 3) = {31, 30, 29};
Physical Line("Inflow", 2) = {44, 48, 50};

Physical Surface("XPML",100) = {78, 68};
Physical Surface("YPML",200) = {82, 64, 72, 74};
Physical Surface("XYPML",300) = {80, 66, 70, 76};
Physical Surface("Interior",9) = {60, 57, 62};

Coherence;