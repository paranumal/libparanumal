r1 = DefineNumber[0.02];
r2 = DefineNumber[0.1];
r3 = DefineNumber[0.2];
r4 = DefineNumber[0.2];
//+
ffxmin = DefineNumber[-0.05];
ffxmax = DefineNumber[0.05];
ffymin = DefineNumber[0.0];
ffymax = DefineNumber[0.5];

//+
fbfacx = DefineNumber[4.0];
fbfacy = DefineNumber[1.5];

fbxmin = DefineNumber[3.0*ffxmin];
fbxmax = DefineNumber[4.0*ffxmax];
fbymin = DefineNumber[1.5*ffymin];
fbymax = DefineNumber[1.5*ffymax];

//+
fdxmin = DefineNumber[-4.0];
fdxmax = DefineNumber[4.0];
fdymin = DefineNumber[0.0];
fdymax = DefineNumber[4.0];

xpmlmin = DefineNumber[-5.0];
xpmlmax = DefineNumber[ 5.0]; 
ypmlmin = DefineNumber[ 0.0];
ypmlmax = DefineNumber[ 5.0];


//+
Point(1) = {ffxmin, ffymin, 0.0, r1};
Point(2) = {ffxmin, ffymax, 0.0, r1};
Point(3) = {ffxmax, ffymax, 0.0, r1};
Point(4) = {ffxmax, ffymin, 0.0, r1};

//+
Point(5) = {fbxmin, fbymin, 0.0, r1};
Point(6) = {fbxmin, fbymax, 0.0, r1};
Point(7) = {fbxmax, fbymax, 0.0, r1};
Point(8) = {fbxmax, fbymin, 0.0, r1};

//+
Point(9)  = {fdxmin, fdymin, 0.0, r3};
Point(10) = {fdxmin, fdymax, 0.0, r3};
Point(11) = {fdxmax, fdymax, 0.0, r3};
Point(12) = {fdxmax, fdymin, 0.0, r3};

//+
Point(17) = {xpmlmin, ypmlmin, 0.0, r3};
Point(18) = {xpmlmax, ypmlmin, 0.0, r3};
Point(19) = {xpmlmin, ypmlmax, 0.0, r3};
Point(20) = {xpmlmax, ypmlmax, 0.0, r3};


Point(21) = {xpmlmin, fdymax, 0.0, r3};
Point(22) = {fdxmin,  ypmlmax, 0.0, r3};

Point(23) = {xpmlmax, fdymax, 0.0, r3};
Point(24) = {fdxmax,  ypmlmax, 0.0, r3};

Point(25) = {-3.0, 0,   0.0, r3};
Point(26) = { 3.0, 0,   0.0, r2};
Point(27) = { 3.0, 2.0, 0.0, r2};
Point(28) = { fdxmax, 2.0, 0.0, r3};

Line(1) = {9, 25};
Line(2) = {25, 5};
Line(3) = {5, 1};
Line(4) = {1, 2};
Line(5) = {2, 3};
Line(6) = {3, 4};
Line(7) = {4, 8};
Line(8) = {8, 26};
Line(9) = {26, 12};
Line(10) = {12, 28};
Line(11) = {28, 27};
Line(12) = {27, 7};
Line(13) = {7, 8};
Line(14) = {7, 6};
Line(15) = {6, 5};
Line(16) = {28, 11};
Line(17) = {11, 10};
Line(18) = {10, 9};
Line(19) = {12, 18};
Line(20) = {18, 23};
Line(21) = {23, 11};
Line(22) = {23, 20};
Line(23) = {20, 24};
Line(24) = {24, 11};
Line(25) = {24, 22};
Line(26) = {22, 10};
Line(27) = {22, 19};
Line(28) = {19, 21};
Line(29) = {21, 10};
Line(30) = {21, 17};
Line(31) = {17, 9};
Line(301) = {27, 26};



Line Loop(302) = {13, 8, -301, 12};
Line Loop(304) = {11, 301, 9, 10};
Plane Surface(305) = {304};
Plane Surface(303) = {302};


//+
Line Loop(32) = {3, 4, 5, 6, 7, -13, 14, 15};
Plane Surface(33) = {32};
//+Line Loop(34) = {13, 8, 9, 10, 11, 12};
//+Plane Surface(35) = {34};
Line Loop(36) = {1, 2, -15, -14, -12, -11, 16, 17, 18};
Plane Surface(37) = {36};
Line Loop(38) = {19, 20, 21, -16, -10};
Plane Surface(39) = {38};
Line Loop(40) = {30, 31, -18, -29};
Plane Surface(41) = {40};
Line Loop(201) = {29, -26, 27, 28};
Plane Surface(202) = {201};
Line Loop(203) = {24, -21, 22, 23};
Plane Surface(204) = {203};
Line Loop(101) = {25, 26, -17, -24};
Plane Surface(102) = {101};


Physical Surface("Interior",9) = {33,305,303,37};
Physical Surface("XPML",100) = {39, 41};
Physical Surface("YPML",200) = {102};
Physical Surface("XYPML",300) = {202,204};

Physical Line("Wall",1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};
Physical Line("Outflow",3) = {19, 20, 22, 23, 25, 27, 28, 30, 31};
//+Physical Line("Slip",4) = {1};
Coherence; 


//+Point{34,35,36} In Surface{38};
//+Field[1] = Attractor;
//+Field[1].FieldX = -1;
//+Field[1].FieldY = -1;
//+Field[1].FieldZ = -1;
//+Field[1].NNodesByEdge = 20;
//+Field[1].NodesList = {6,7,8,26,34};
//+Field[2] = Threshold;
//+Field[2].DistMax = 2;
//+Field[2].DistMin = 0.5;
//+Field[2].IField = 1;
//+Field[2].LcMax = 0.2;
//+Field[2].LcMin = 0.01;
//+Field[2].Sigmoid = 0;
//+Field[2].StopAtDistMax = 0;
//+Field[3] = Min;
//+Field[3].FieldsList = {2};
//+Background Field = 3;
