// K=18,832 [quads]

lc = 1.00;

// x-bounds
xmin  = -10.0;
pmlx1 =  -8.0;
// pmlx2 =  18.0;
// xmax  =  20.0;
pmlx2 =  13.0;
xmax  =  15.0;


// y-bounds
ymin  =  0.0;
pmly1 =  0.0;
pmly2 =  6.0;
ymax  =  8.0;

// z-bounds
zall = 0.0;

fx1 = 0.00;         // front of fence
fx2 = fx1 + 0.10;   // width of fence
fy1 = 0.00;         // base of fence
fy2 = 1.00;         // top of fence


Point(1) = { xmin,   ymin, zall, lc*1.0};
Point(2) = {pmlx1,   ymin, zall, lc*1.0};
//---------------------------------------
Point(3) = { fx1,  fy1, zall, 0.050};  // base front
Point(4) = { fx1,  fy2, zall, 0.025};  // top  front
Point(5) = { fx2,  fy2, zall, 0.025};  // top  rear
Point(6) = { fx2,  fy1, zall, 0.050};  // base rear
//---------------------------------------
Point(7) = {pmlx2,  ymin,  zall, lc*1.0};
Point(8) = {xmax,   ymin,  zall, lc*1.0};

Point( 9)= {xmax,   pmly2, zall, lc*1.0};
Point(10)= {xmax,   ymax,  zall, lc*1.0};
Point(11)= {pmlx2,  ymax,  zall, lc*1.0};
Point(12)= {pmlx1,  ymax,  zall, lc*1.0};
Point(13)= {xmin,   ymax,  zall, lc*1.0};
//---------------------------------------
Point(14)= {pmlx2, pmly2, zall, lc*1.0};
Point(15)= {pmlx1, pmly2, zall, lc*1.0};
Point(16)= {xmin,  pmly2, zall, lc*1.0};
//---------------------------------------
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,9};

Line(11)= { 9,10};
Line(12)= {10,11};
Line(13)= {11,12};
Line(14)= {12,13};
Line(15)= {13,16};
Line(16)= {16, 1};
Line(17)= {15, 2};
Line(18)= {15,12};
Line(19)= { 7,14};
Line(20)= {14,11};
Line(21)= { 9,14};
Line(22)= {14,15};
Line(23)= {15,16};


//---------------------------------------
// main fluid domain
//---------------------------------------
LLDomain=newll; Line Loop(LLDomain) = {2,3,4,5,6,19,22,17};
PSDomain=news;  Plane Surface(PSDomain) = {LLDomain};

//---------------------------------------
// 5 PML regions
//---------------------------------------
LLpmlA=newll;     Line Loop(LLpmlA) = {1,-17,23,16};
PSpmlA=news;  Plane Surface(PSpmlA) = {LLpmlA};

LLpmlB=newll;     Line Loop(LLpmlB) = {-23,18,14,15};
PSpmlB=news;  Plane Surface(PSpmlB) = {LLpmlB};

LLpmlC=newll;     Line Loop(LLpmlC) = {-22,20,13,-18};
PSpmlC=news;  Plane Surface(PSpmlC) = {LLpmlC};

LLpmlD=newll;     Line Loop(LLpmlD) = {-21,11,12,-20};
PSpmlD=news;  Plane Surface(PSpmlD) = {LLpmlD};

LLpmlE=newll;     Line Loop(LLpmlE) = {7,8,21,-19};
PSpmlE=news;  Plane Surface(PSpmlE) = {LLpmlE};

Physical Surface("Domian",9) = {PSDomain};
Physical Surface("PMLX",100) = {PSpmlA,PSpmlE};
Physical Surface("PMLY",200) = {PSpmlC};
Physical Surface("PMLXY",300) = {PSpmlB,PSpmlD};
Physical    Line("Wall", 1) = {2,3,4,5,6};
Physical    Line("Inflow",2) = {15,16};
Physical    Line("Outflow",3) = {8,11};
Physical    Line("Slip",4) = {1,7,12,13,14};

//---------------------------------------------------------
// use Fields to tune wake density
//---------------------------------------------------------

Field[1] = Box;
Field[1].VIn  =  0.035;
Field[1].VOut =  0.50;
Field[1].XMin = -0.8;
Field[1].XMax =  2.5;
Field[1].YMin =  0.0;
Field[1].YMax =  1.8;

Field[2] = Box;
Field[2].VIn  =  0.05;
Field[2].VOut =  0.50;
Field[2].XMin = -2.2;
Field[2].XMax =  5.0;
Field[2].YMin =  0.0;
Field[2].YMax =  2.7;

Field[3] = Box;
Field[3].VIn  =  0.05;
Field[3].VOut =  0.50;
Field[3].XMin = -3.3;
Field[3].XMax =  0.0;
Field[3].YMin =  0.0;
Field[3].YMax =  1.0;

Field[4] = Box;
Field[4].VIn  =  0.10;
Field[4].VOut =  0.50;
Field[4].XMin = -4.5;
Field[4].XMax =  8.0;
Field[4].YMin =  0.0;
Field[4].YMax =  3.0;

Field[5] = Box;
Field[5].VIn  =  0.15;
Field[5].VOut =  0.50;
Field[5].XMin = -6.0;
Field[5].XMax = 12.0;
Field[5].YMin =  0.0;
Field[5].YMax =  4.0;

Field[6] = Box;
Field[6].VIn  =  0.02;
Field[6].VOut =  0.50;
Field[6].XMin = -0.2;
Field[6].XMax =  0.5;
Field[6].YMin =  0.8;
Field[6].YMax =  1.2;


Field[7] = Min;
Field[7].FieldsList = {1,2,3,4,5,6};
Background Field = 7;
