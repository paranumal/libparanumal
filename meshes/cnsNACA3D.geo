Geometry.OldNewReg=0;
SetFactory("OpenCASCADE");
// Mesh.LcIntegrationPrecision = 1.e-2;

lc   = DefineNumber[0.05];
lc2  = DefineNumber[0.25];
lc3  = DefineNumber[0.5];

Point(1) =  {1.000000e+00,0.000000e+00,0.000000e+00, lc};
Point(2) =  {9.997533e-01,0.000000e+00,-3.498543e-05,lc};
Point(3) =  {9.990134e-01,0.000000e+00,-1.398841e-04,lc};
Point(4) =  {9.977810e-01,0.000000e+00,-3.143904e-04,lc};
Point(5) =  {9.960575e-01,0.000000e+00,-5.579769e-04,lc};
Point(6) =  {9.938442e-01,0.000000e+00,-8.699747e-04,lc};
Point(7) =  {9.911436e-01,0.000000e+00,-1.249551e-03,lc};
Point(8) =  {9.879584e-01,0.000000e+00,-1.695579e-03,lc};
Point(9) =  {9.842916e-01,0.000000e+00,-2.206860e-03,lc};
Point(10) = {9.801469e-01,0.000000e+00,-2.781989e-03,lc};
Point(11) = {9.755284e-01,0.000000e+00,-3.419365e-03,lc};
Point(12) = {9.704404e-01,0.000000e+00,-4.117359e-03,lc};
Point(13) = {9.648883e-01,0.000000e+00,-4.874101e-03,lc};
Point(14) = {9.588774e-01,0.000000e+00,-5.687566e-03,lc};
Point(15) = {9.524136e-01,0.000000e+00,-6.555737e-03,lc};
Point(16) = {9.455034e-01,0.000000e+00,-7.476377e-03,lc};
Point(17) = {9.381535e-01,0.000000e+00,-8.447210e-03,lc};
Point(18) = {9.303712e-01,0.000000e+00,-9.465891e-03,lc};
Point(19) = {9.221641e-01,0.000000e+00,-1.052998e-02,lc};
Point(20) = {9.135405e-01,0.000000e+00,-1.163695e-02,lc};
Point(21) = {9.045087e-01,0.000000e+00,-1.278429e-02,lc};
Point(22) = {8.950777e-01,0.000000e+00,-1.396934e-02,lc};
Point(23) = {8.852569e-01,0.000000e+00,-1.518951e-02,lc};
Point(24) = {8.750558e-01,0.000000e+00,-1.644214e-02,lc};
Point(25) = {8.644845e-01,0.000000e+00,-1.772453e-02,lc};
Point(26) = {8.535537e-01,0.000000e+00,-1.903398e-02,lc};
Point(27) = {8.422739e-01,0.000000e+00,-2.036772e-02,lc};
Point(28) = {8.306563e-01,0.000000e+00,-2.172309e-02,lc};
Point(29) = {8.187124e-01,0.000000e+00,-2.309725e-02,lc};
Point(30) = {8.064539e-01,0.000000e+00,-2.448751e-02,lc};
Point(31) = {7.938930e-01,0.000000e+00,-2.589105e-02,lc};
Point(32) = {7.810421e-01,0.000000e+00,-2.730503e-02,lc};
Point(33) = {7.679139e-01,0.000000e+00,-2.872668e-02,lc};
Point(34) = {7.545212e-01,0.000000e+00,-3.015313e-02,lc};
Point(35) = {7.408773e-01,0.000000e+00,-3.158154e-02,lc};
Point(36) = {7.269957e-01,0.000000e+00,-3.300894e-02,lc};
Point(37) = {7.128901e-01,0.000000e+00,-3.443245e-02,lc};
Point(38) = {6.985745e-01,0.000000e+00,-3.584905e-02,lc};
Point(39) = {6.840628e-01,0.000000e+00,-3.725576e-02,lc};
Point(40) = {6.693696e-01,0.000000e+00,-3.864942e-02,lc};
Point(41) = {6.545091e-01,0.000000e+00,-4.002701e-02,lc};
Point(42) = {6.394961e-01,0.000000e+00,-4.138529e-02,lc};
Point(43) = {6.243456e-01,0.000000e+00,-4.272101e-02,lc};
Point(44) = {6.090723e-01,0.000000e+00,-4.403092e-02,lc};
Point(45) = {5.936913e-01,0.000000e+00,-4.531165e-02,lc};
Point(46) = {5.782179e-01,0.000000e+00,-4.655984e-02,lc};
Point(47) = {5.626673e-01,0.000000e+00,-4.777199e-02,lc};
Point(48) = {5.470549e-01,0.000000e+00,-4.894463e-02,lc};
Point(49) = {5.313960e-01,0.000000e+00,-5.007425e-02,lc};
Point(50) = {5.157061e-01,0.000000e+00,-5.115728e-02,lc};
Point(51) = {5.000008e-01,0.000000e+00,-5.219014e-02,lc};
Point(52) = {4.842954e-01,0.000000e+00,-5.316926e-02,lc};
Point(53) = {4.686055e-01,0.000000e+00,-5.409108e-02,lc};
Point(54) = {4.529467e-01,0.000000e+00,-5.495201e-02,lc};
Point(55) = {4.373342e-01,0.000000e+00,-5.574857e-02,lc};
Point(56) = {4.217836e-01,0.000000e+00,-5.647729e-02,lc};
Point(57) = {4.063102e-01,0.000000e+00,-5.713477e-02,lc};
Point(58) = {3.909292e-01,0.000000e+00,-5.771770e-02,lc};
Point(59) = {3.756559e-01,0.000000e+00,-5.822293e-02,lc};
Point(60) = {3.605053e-01,0.000000e+00,-5.864737e-02,lc};
Point(61) = {3.454924e-01,0.000000e+00,-5.898812e-02,lc};
Point(62) = {3.306319e-01,0.000000e+00,-5.924247e-02,lc};
Point(63) = {3.159386e-01,0.000000e+00,-5.940786e-02,lc};
Point(64) = {3.014269e-01,0.000000e+00,-5.948193e-02,lc};
Point(65) = {2.871112e-01,0.000000e+00,-5.946260e-02,lc};
Point(66) = {2.730056e-01,0.000000e+00,-5.934800e-02,lc};
Point(67) = {2.591240e-01,0.000000e+00,-5.913650e-02,lc};
Point(68) = {2.454802e-01,0.000000e+00,-5.882679e-02,lc};
Point(69) = {2.320875e-01,0.000000e+00,-5.841779e-02,lc};
Point(70) = {2.189592e-01,0.000000e+00,-5.790876e-02,lc};
Point(71) = {2.061082e-01,0.000000e+00,-5.729925e-02,lc};
Point(72) = {1.935473e-01,0.000000e+00,-5.658907e-02,lc};
Point(73) = {1.812888e-01,0.000000e+00,-5.577839e-02,lc};
Point(74) = {1.693449e-01,0.000000e+00,-5.486767e-02,lc};
Point(75) = {1.577273e-01,0.000000e+00,-5.385765e-02,lc};
Point(76) = {1.464474e-01,0.000000e+00,-5.274938e-02,lc};
Point(77) = {1.355165e-01,0.000000e+00,-5.154420e-02,lc};
Point(78) = {1.249452e-01,0.000000e+00,-5.024372e-02,lc};
Point(79) = {1.147441e-01,0.000000e+00,-4.884978e-02,lc};
Point(80) = {1.049232e-01,0.000000e+00,-4.736451e-02,lc};
Point(81) = {9.549212e-02,0.000000e+00,-4.579021e-02,lc};
Point(82) = {8.646032e-02,0.000000e+00,-4.412942e-02,lc};
Point(83) = {7.783660e-02,0.000000e+00,-4.238483e-02,lc};
Point(84) = {6.962952e-02,0.000000e+00,-4.055926e-02,lc};
Point(85) = {6.184718e-02,0.000000e+00,-3.865567e-02,lc};
Point(86) = {5.449721e-02,0.000000e+00,-3.667711e-02,lc};
Point(87) = {4.758692e-02,0.000000e+00,-3.462668e-02,lc};
Point(88) = {4.112309e-02,0.000000e+00,-3.250752e-02,lc};
Point(89) = {3.511214e-02,0.000000e+00,-3.032277e-02,lc};
Point(90) = {2.955997e-02,0.000000e+00,-2.807550e-02,lc};
Point(91) = {2.447206e-02,0.000000e+00,-2.576878e-02,lc};
Point(92) = {1.985344e-02,0.000000e+00,-2.340553e-02,lc};
Point(93) = {1.570869e-02,0.000000e+00,-2.098859e-02,lc};
Point(94) = {1.204184e-02,0.000000e+00,-1.852062e-02,lc};
Point(95) = {8.856565e-03,0.000000e+00,-1.600414e-02,lc};
Point(96) = {6.155997e-03,0.000000e+00,-1.344148e-02,lc};
Point(97) = {3.942788e-03,0.000000e+00,-1.083471e-02,lc};
Point(98) = {2.219111e-03,0.000000e+00,-8.185687e-03,lc};
Point(99) = {9.866953e-04,0.000000e+00,-5.496060e-03,lc};
Point(100) = {2.467632e-04,0.000000e+00,-2.767267e-03,lc};

Point(101) = {0.000000e+00,0.000000e+00,1.911503e-39,lc};
Point(102) = {2.467632e-04,0.000000e+00,2.767267e-03,lc};
Point(103) = {9.866953e-04,0.000000e+00,5.496060e-03,lc};
Point(104) = {2.219111e-03,0.000000e+00,8.185687e-03,lc};
Point(105) = {3.942788e-03,0.000000e+00,1.083471e-02,lc};
Point(106) = {6.155997e-03,0.000000e+00,1.344148e-02,lc};
Point(107) = {8.856565e-03,0.000000e+00,1.600414e-02,lc};
Point(108) = {1.204184e-02,0.000000e+00,1.852062e-02,lc};
Point(109) = {1.570869e-02,0.000000e+00,2.098859e-02,lc};
Point(110) = {1.985344e-02,0.000000e+00,2.340553e-02,lc};
Point(111) = {2.447206e-02,0.000000e+00,2.576878e-02,lc};
Point(112) = {2.955997e-02,0.000000e+00,2.807550e-02,lc};
Point(113) = {3.511214e-02,0.000000e+00,3.032277e-02,lc};
Point(114) = {4.112309e-02,0.000000e+00,3.250752e-02,lc};
Point(115) = {4.758692e-02,0.000000e+00,3.462668e-02,lc};
Point(116) = {5.449721e-02,0.000000e+00,3.667711e-02,lc};
Point(117) = {6.184718e-02,0.000000e+00,3.865567e-02,lc};
Point(118) = {6.962952e-02,0.000000e+00,4.055926e-02,lc};
Point(119) = {7.783660e-02,0.000000e+00,4.238483e-02,lc};
Point(120) = {8.646032e-02,0.000000e+00,4.412942e-02,lc};
Point(121) = {9.549212e-02,0.000000e+00,4.579021e-02,lc};
Point(122) = {1.049232e-01,0.000000e+00,4.736451e-02,lc};
Point(123) = {1.147441e-01,0.000000e+00,4.884978e-02,lc};
Point(124) = {1.249452e-01,0.000000e+00,5.024372e-02,lc};
Point(125) = {1.355165e-01,0.000000e+00,5.154420e-02,lc};
Point(126) = {1.464474e-01,0.000000e+00,5.274938e-02,lc};
Point(127) = {1.577273e-01,0.000000e+00,5.385765e-02,lc};
Point(128) = {1.693449e-01,0.000000e+00,5.486767e-02,lc};
Point(129) = {1.812888e-01,0.000000e+00,5.577839e-02,lc};
Point(130) = {1.935473e-01,0.000000e+00,5.658907e-02,lc};
Point(131) = {2.061082e-01,0.000000e+00,5.729925e-02,lc};
Point(132) = {2.189592e-01,0.000000e+00,5.790876e-02,lc};
Point(133) = {2.320875e-01,0.000000e+00,5.841779e-02,lc};
Point(134) = {2.454802e-01,0.000000e+00,5.882679e-02,lc};
Point(135) = {2.591240e-01,0.000000e+00,5.913650e-02,lc};
Point(136) = {2.730056e-01,0.000000e+00,5.934800e-02,lc};
Point(137) = {2.871112e-01,0.000000e+00,5.946260e-02,lc};
Point(138) = {3.014269e-01,0.000000e+00,5.948193e-02,lc};
Point(139) = {3.159386e-01,0.000000e+00,5.940786e-02,lc};
Point(140) = {3.306319e-01,0.000000e+00,5.924247e-02,lc};
Point(141) = {3.454924e-01,0.000000e+00,5.898812e-02,lc};
Point(142) = {3.605053e-01,0.000000e+00,5.864737e-02,lc};
Point(143) = {3.756559e-01,0.000000e+00,5.822293e-02,lc};
Point(144) = {3.909292e-01,0.000000e+00,5.771770e-02,lc};
Point(145) = {4.063102e-01,0.000000e+00,5.713477e-02,lc};
Point(146) = {4.217836e-01,0.000000e+00,5.647729e-02,lc};
Point(147) = {4.373342e-01,0.000000e+00,5.574857e-02,lc};
Point(148) = {4.529467e-01,0.000000e+00,5.495201e-02,lc};
Point(149) = {4.686055e-01,0.000000e+00,5.409108e-02,lc};
Point(150) = {4.842954e-01,0.000000e+00,5.316926e-02,lc};
Point(151) = {5.000008e-01,0.000000e+00,5.219014e-02,lc};
Point(152) = {5.157061e-01,0.000000e+00,5.115728e-02,lc};
Point(153) = {5.313960e-01,0.000000e+00,5.007425e-02,lc};
Point(154) = {5.470549e-01,0.000000e+00,4.894463e-02,lc};
Point(155) = {5.626673e-01,0.000000e+00,4.777199e-02,lc};
Point(156) = {5.782179e-01,0.000000e+00,4.655984e-02,lc};
Point(157) = {5.936913e-01,0.000000e+00,4.531165e-02,lc};
Point(158) = {6.090723e-01,0.000000e+00,4.403092e-02,lc};
Point(159) = {6.243456e-01,0.000000e+00,4.272101e-02,lc};
Point(160) = {6.394961e-01,0.000000e+00,4.138529e-02,lc};
Point(161) = {6.545091e-01,0.000000e+00,4.002701e-02,lc};
Point(162) = {6.693696e-01,0.000000e+00,3.864942e-02,lc};
Point(163) = {6.840628e-01,0.000000e+00,3.725576e-02,lc};
Point(164) = {6.985745e-01,0.000000e+00,3.584905e-02,lc};
Point(165) = {7.128901e-01,0.000000e+00,3.443245e-02,lc};
Point(166) = {7.269957e-01,0.000000e+00,3.300894e-02,lc};
Point(167) = {7.408773e-01,0.000000e+00,3.158154e-02,lc};
Point(168) = {7.545212e-01,0.000000e+00,3.015313e-02,lc};
Point(169) = {7.679139e-01,0.000000e+00,2.872668e-02,lc};
Point(170) = {7.810421e-01,0.000000e+00,2.730503e-02,lc};
Point(171) = {7.938930e-01,0.000000e+00,2.589105e-02,lc};
Point(172) = {8.064539e-01,0.000000e+00,2.448751e-02,lc};
Point(173) = {8.187124e-01,0.000000e+00,2.309725e-02,lc};
Point(174) = {8.306563e-01,0.000000e+00,2.172309e-02,lc};
Point(175) = {8.422739e-01,0.000000e+00,2.036772e-02,lc};
Point(176) = {8.535537e-01,0.000000e+00,1.903398e-02,lc};
Point(177) = {8.644845e-01,0.000000e+00,1.772453e-02,lc};
Point(178) = {8.750558e-01,0.000000e+00,1.644214e-02,lc};
Point(179) = {8.852569e-01,0.000000e+00,1.518951e-02,lc};
Point(180) = {8.950777e-01,0.000000e+00,1.396934e-02,lc};
Point(181) = {9.045087e-01,0.000000e+00,1.278429e-02,lc};
Point(182) = {9.135405e-01,0.000000e+00,1.163695e-02,lc};
Point(183) = {9.221641e-01,0.000000e+00,1.052998e-02,lc};
Point(184) = {9.303712e-01,0.000000e+00,9.465891e-03,lc};
Point(185) = {9.381535e-01,0.000000e+00,8.447210e-03,lc};
Point(186) = {9.455034e-01,0.000000e+00,7.476377e-03,lc};
Point(187) = {9.524136e-01,0.000000e+00,6.555737e-03,lc};
Point(188) = {9.588774e-01,0.000000e+00,5.687566e-03,lc};
Point(189) = {9.648883e-01,0.000000e+00,4.874101e-03,lc};
Point(190) = {9.704404e-01,0.000000e+00,4.117359e-03,lc};
Point(191) = {9.755284e-01,0.000000e+00,3.419365e-03,lc};
Point(192) = {9.801469e-01,0.000000e+00,2.781989e-03,lc};
Point(193) = {9.842916e-01,0.000000e+00,2.206860e-03,lc};
Point(194) = {9.879584e-01,0.000000e+00,1.695579e-03,lc};
Point(195) = {9.911436e-01,0.000000e+00,1.249551e-03,lc};
Point(196) = {9.938442e-01,0.000000e+00,8.699747e-04,lc};
Point(197) = {9.960575e-01,0.000000e+00,5.579769e-04,lc};
Point(198) = {9.977810e-01,0.000000e+00,3.143904e-04,lc};
Point(199) = {9.990134e-01,0.000000e+00,1.398841e-04,lc};
Point(200) = {9.997533e-01,0.000000e+00,3.498543e-05,lc};

Spline(1) = { 1 ... 50};
Spline(2) = { 50 ... 101};
Spline(3) = { 101 ... 151};
Spline(4) = { 151 ... 200,1};


// move to x-y plane
Rotate { {1,0,0},{0,0,0},Pi/2 } { Line{1,2,3,4}; }

// move center to 0.0.0 point

L = 1.0; 
Translate {-0.5,0,-L/2} { Line{1,2,3,4};}

// Create Wing
Curve Loop(1) = {1,2,3,4}; 
Plane Surface(1) = {1};
Extrude {0, 0, L} { Surface{1};}

// n = #out[];
// Printf("Extrude (line) has returned %g elements", n);
// Printf("Wing Extrude Returned Volume Id: %d ", out[0]);

// Printf("next point ID = %g", newp);
// Printf("next point ID = %g; next line ID = %g; next surface ID = %g", newp,  newl, news);

xc   = DefineNumber[0];
yc   = DefineNumber[0];
zc   = DefineNumber[-L/2];

xi   = DefineNumber[0];
yi   = DefineNumber[1.0];
zi   = DefineNumber[-L/2];

yomin   = DefineNumber[-3.0];
yomax   = -yomin; 
zomin   = -L/2; 

xomin   = DefineNumber[3.0];
xomax   = DefineNumber[4.0];


// External box
p1 = newp; Point(p1) = { xc, yc,  zomin, lc3}; 
p2 = newp; Point(p2) = { xc, yomin,  zomin, lc3}; 
p3 = newp; Point(p3) = { xomax, yomin,  zomin, lc3}; 
p4 = newp; Point(p4) = { xomax, yomax,  zomin, lc3}; 
p5 = newp; Point(p5) = { xc,   yomax,  zomin, lc3}; 
p6 = newp; Point(p6) = { yomin,  yc,  zomin, lc3}; 

// Internal region
p7 = newp; Point(p7) = { xi, -yi,  zomin, lc2}; 
p8 = newp; Point(p8) = { xi,  yi,  zomin, lc2}; 
p9 = newp; Point(p9) = { -yi,  0.0,  zomin, lc2}; 

p10 = newp;  Point(p10) = { xomax, -(yi+0.5),  zomin, lc2}; 
p11 = newp;  Point(p11) = { xomax,  (yi+0.5),  zomin, lc2}; 

l1 = newc; Line(l1) = {p2,  p3};
l2 = newc; Line(l2) = {p3,  p10};
l3 = newc; Line(l3) = {p10, p11};
l4 = newc; Line(l4) = {p11, p4}; 
l5 = newc; Line(l5) = {p4,  p5}; 

c1 = newc; Circle(c1) = {p5, p1, p6};
c2 = newc; Circle(c2) = {p6, p1, p2};

l6 = newc; Line(l6)  = {p7,p10}; 
l7 = newc; Line(l7)  = {p8,p11}; 

c3 = newc; Circle(c3) = {p8, p1, p9};
c4 = newc; Circle(c4) = {p9, p1, p7};

//+
Curve Loop(7) = {22, 23, 20, 15, -21};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {18, 19, 13, 14, -20, -23, -22, 21, 16, 17};
//+
Plane Surface(8) = {8};
Extrude {0, 0, 3*L} { Surface{7,8};}

BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; }

Physical Surface("Farfield", 20) = {15, 16, 17, 20, 19, 24, 18, 21, 23};
Physical Surface("Wall", 11) = {25, 26, 29, 27, 28};
Physical Surface("Symmetry", 13) = {8, 22};
Physical Volume("Domain", 9) = {3, 2};


Coherence;


Mesh.MeshSizeFromCurvature = 20;
