rf = DefineNumber[0.25];
rd = DefineNumber[0.75];
// trailing side resolution
rf2 = DefineNumber[0.5];


ffxmin = DefineNumber[-0.5];
ffxmax = DefineNumber[0.5];
ffymin = DefineNumber[0.0];
ffymax = DefineNumber[2.0];
ffzmin = DefineNumber[-0.5];
ffzmax = DefineNumber[0.5];
//+
fbfac2 = DefineNumber[1.5];
fbfac1 = DefineNumber[3];

fbxmin = DefineNumber[fbfac1*ffxmin];
fbxmax = DefineNumber[fbfac1*ffxmax];
fbymin = DefineNumber[fbfac2*ffymin];
fbymax = DefineNumber[fbfac2*ffymax];
fbzmin = DefineNumber[fbfac1*ffzmin];
fbzmax = DefineNumber[fbfac1*ffzmax];


fdxmin = DefineNumber[-6];
fdxmax = DefineNumber[12];
fdymin = DefineNumber[0.0];
fdymax = DefineNumber[6.0];
fdzmin = DefineNumber[-4.5];
fdzmax = DefineNumber[4.5];



//+
Point(1) = {ffxmin, ffymin, ffzmin, rf};
Point(2) = {ffxmin, ffymax, ffzmin, rf};
Point(3) = {ffxmax, ffymax, ffzmin, rf};
Point(4) = {ffxmax, ffymin, ffzmin, rf};
//+
Point(5) = {ffxmin, ffymin, ffzmax, rf};
Point(6) = {ffxmin, ffymax, ffzmax, rf};
Point(7) = {ffxmax, ffymax, ffzmax, rf};
Point(8) = {ffxmax, ffymin, ffzmax, rf};

Line(1) = {5, 8};
Line(2) = {8, 7};
Line(3) = {7, 6};
Line(4) = {6, 5};
Line(5) = {6, 2};
Line(6) = {2, 1};
Line(7) = {1, 5};
Line(8) = {2, 3};
Line(9) = {3, 4};
Line(10) = {4, 1};
Line(11) = {3, 7};
Line(12) = {4, 8};
Line Loop(13) = {8, 11, 3, 5};
Plane Surface(14) = {13};
Line Loop(15) = {8, 9, 10, -6};
Plane Surface(16) = {15};
Line Loop(17) = {9, 12, 2, -11};
Plane Surface(18) = {17};
Line Loop(19) = {3, 4, 1, 2};
Plane Surface(20) = {19};
Line Loop(21) = {7, -4, 5, 6};
Plane Surface(22) = {21};
//+

//+
Point(9)  = {fbxmin, fbymin, fbzmin, rf};
Point(10) = {fbxmin, fbymax, fbzmin, rf};
Point(11) = {fbxmax, fbymax, fbzmin, rf};
Point(12) = {fbxmax, fbymin, fbzmin, rf};
//+
Point(13) = {fbxmin, fbymin, fbzmax, rf};
Point(14) = {fbxmin, fbymax, fbzmax, rf};
Point(15) = {fbxmax, fbymax, fbzmax, rf};
Point(16) = {fbxmax, fbymin, fbzmax, rf};

Line(23) = {13, 16};
Line(24) = {16, 12};
Line(25) = {12, 9};
Line(26) = {9, 13};
Line(27) = {16, 15};
Line(28) = {15, 14};
Line(29) = {14, 13};
Line(30) = {14, 10};
Line(31) = {10, 9};
Line(32) = {10, 11};
Line(33) = {11, 12};
Line(34) = {11, 15};
Line Loop(35) = {26, 23, 24, 25};
Line Loop(36) = {1, -12, 10, 7};
Plane Surface(37) = {35, 36};
Line Loop(38) = {25, -31, 32, 33};
Plane Surface(39) = {38};
Line Loop(40) = {33, -24, 27, -34};
Plane Surface(41) = {40};
Line Loop(42) = {32, 34, 28, 30};
Plane Surface(43) = {42};
Line Loop(44) = {27, 28, 29, 23};
Plane Surface(45) = {44};
Line Loop(46) = {30, 31, 26, -29};
Plane Surface(47) = {46};
Surface Loop(48) = {37, 47, 43, 39, 41, 45, 22, 20, 14, 16, 18};
Volume(49) = {48};
//+






Point(17) = {fdxmin, fbymin, fbzmin, rd};
Point(18) = {fdxmin, fbymax, fbzmin, rd};
Point(19) = {fdxmax, fbymax, fbzmin, rf2};
Point(20) = {fdxmax, fbymin, fbzmin, rf2};

Point(21) = {fdxmin, fbymin, fbzmax, rd};
Point(22) = {fdxmin, fbymax, fbzmax, rd};
Point(23) = {fdxmax, fbymax, fbzmax, rf2};
Point(24) = {fdxmax, fbymin, fbzmax, rf2};
Line(50) = {14, 22};
Line(51) = {22, 21};
Line(52) = {21, 13};
Line(53) = {22, 18};
Line(54) = {18, 17};
Line(55) = {17, 21};
Line(56) = {17, 9};
Line(57) = {10, 18};
Line(58) = {19, 11};
Line(59) = {19, 20};
Line(60) = {20, 12};
Line(61) = {15, 23};
Line(62) = {19, 23};
Line(63) = {23, 24};
Line(64) = {24, 20};
Line(65) = {24, 16};
Line Loop(66) = {50, 53, -57, -30};
Plane Surface(67) = {66};
Line Loop(68) = {51, 52, -29, 50};
Plane Surface(69) = {68};
Line Loop(70) = {51, -55, -54, -53};
Plane Surface(71) = {70};
Line Loop(72) = {57, 54, 56, -31};
Plane Surface(73) = {72};
Line Loop(74) = {55, 52, -26, -56};
Plane Surface(75) = {74};
Surface Loop(76) = {67, 69, 71, 75, 73, 39, 37, 45, 41, 43, 22, 20, 14, 16, 18};
Surface Loop(77) = {67, 69, 71, 75, 73, 47};
Volume(78) = {77};
Line Loop(79) = {61, -62, 58, 34};
Plane Surface(80) = {79};
Line Loop(81) = {27, 61, 63, 65};
Plane Surface(82) = {81};
Line Loop(83) = {24, -60, -64, 65};
Plane Surface(84) = {83};
Line Loop(85) = {62, 63, 64, -59};
Plane Surface(86) = {85};
Line Loop(87) = {58, 33, -60, -59};
Plane Surface(88) = {87};
Surface Loop(89) = {80, 82, 86, 84, 88, 41};
Volume(90) = {89};
//+

Point(25) = {fdxmin, fbymin, fdzmin, rd};
Point(26) = {fdxmin, fbymax, fdzmin, rd};

Point(27) = {fdxmax, fbymax, fdzmin, rd};
Point(28) = {fdxmax, fbymin, fdzmin, rd};

Point(29) = {fbxmax, fbymax, fdzmin, rd};
Point(30) = {fbxmax, fbymin, fdzmin, rd};

Point(31) = {fbxmin, fbymax, fdzmin, rd};
Point(32) = {fbxmin, fbymin, fdzmin, rd};
Line(91) = {20, 28};
Line(92) = {28, 27};
Line(93) = {27, 19};
Line(94) = {27, 29};
Line(95) = {29, 30};
Line(96) = {30, 28};
Line(97) = {30, 12};
Line(98) = {29, 11};
Line(99) = {29, 31};
Line(100) = {31, 32};
Line(101) = {32, 30};
Line(102) = {32, 9};
Line(103) = {10, 31};
Line(104) = {31, 26};
Line(105) = {26, 25};
Line(106) = {25, 32};
Line(107) = {25, 17};
Line(108) = {18, 26};
Line Loop(109) = {93, 59, 91, 92};
Plane Surface(110) = {109};
Line Loop(111) = {92, 94, 95, 96};
Plane Surface(112) = {111};
Line Loop(113) = {96, -91, 60, -97};
Plane Surface(114) = {113};
Line Loop(115) = {95, 97, -33, -98};
Plane Surface(116) = {115};
Line Loop(117) = {58, -98, -94, 93};
Plane Surface(118) = {117};
Surface Loop(119) = {118, 116, 112, 110, 114, 88};
Volume(120) = {119};
Line Loop(121) = {95, -101, -100, -99};
Plane Surface(122) = {121};
Line Loop(123) = {101, 97, 25, -102};
Plane Surface(124) = {123};
Line Loop(125) = {99, -103, 32, -98};
Plane Surface(126) = {125};
Line Loop(127) = {103, 100, 102, -31};
Plane Surface(128) = {127};
Surface Loop(129) = {126, 122, 124, 128, 39, 116};
Volume(130) = {129};
Line Loop(131) = {106, 102, -56, -107};
Plane Surface(132) = {131};
Line Loop(133) = {106, -100, 104, 105};
Plane Surface(134) = {133};
Line Loop(135) = {103, 104, -108, -57};
Plane Surface(136) = {135};
Line Loop(137) = {54, -107, -105, -108};
Plane Surface(138) = {137};
Surface Loop(139) = {136, 134, 132, 138, 128, 73};
Volume(140) = {139};

//+
Point(33) = {fdxmin, fbymin, fdzmax, rd};
Point(34) = {fdxmin, fbymax, fdzmax, rd};

Point(35) = {fdxmax, fbymax, fdzmax, rd};
Point(36) = {fdxmax, fbymin, fdzmax, rd};

Point(37) = {fbxmax, fbymax, fdzmax, rd};
Point(38) = {fbxmax, fbymin, fdzmax, rd};

Point(39) = {fbxmin, fbymax, fdzmax, rd};
Point(40) = {fbxmin, fbymin, fdzmax, rd};

Line(141) = {23, 35};
Line(142) = {35, 36};
Line(143) = {36, 24};
Line(144) = {35, 37};
Line(145) = {37, 38};
Line(146) = {38, 36};
Line(147) = {38, 16};
Line(148) = {15, 37};
Line(149) = {38, 40};
Line(150) = {40, 39};
Line(151) = {39, 37};
Line(152) = {39, 14};
Line(153) = {13, 40};
Line(154) = {40, 33};
Line(155) = {33, 34};
Line(156) = {34, 39};
Line(157) = {34, 22};
Line(158) = {21, 33};
Line Loop(159) = {141, 142, 143, -63};
Plane Surface(160) = {159};
Line Loop(161) = {142, -146, -145, -144};
Plane Surface(162) = {161};
Line Loop(163) = {144, -148, 61, 141};
Plane Surface(164) = {163};
Line Loop(165) = {143, 65, -147, 146};
Plane Surface(166) = {165};
Line Loop(167) = {145, 147, 27, 148};
Plane Surface(168) = {167};
Surface Loop(169) = {164, 162, 160, 166, 168, 82};
Volume(170) = {169};
Line Loop(171) = {145, 149, 150, 151};
Plane Surface(172) = {171};
Line Loop(173) = {147, -23, 153, -149};
Plane Surface(174) = {173};
Line Loop(175) = {148, -151, 152, -28};
Plane Surface(176) = {175};
Plane Surface(177) = {173};
Line Loop(178) = {152, 29, 153, 150};
Plane Surface(179) = {178};
Surface Loop(180) = {176, 172, 179, 168, 45, 174};
Volume(181) = {180};
Line Loop(182) = {156, -150, 154, 155};
Plane Surface(183) = {182};
Line Loop(184) = {156, 152, 50, -157};
Plane Surface(185) = {184};
Line Loop(186) = {153, 154, -158, 52};
Plane Surface(187) = {186};
Line Loop(188) = {157, 51, 158, 155};
Plane Surface(189) = {188};
Surface Loop(190) = {187, 183, 185, 189, 69, 179};
Volume(191) = {190};

//+

//+
Point(41) = {fdxmin, fdymax, fdzmax, rd};
Point(42) = {fdxmax, fdymax, fdzmax, rd};
Point(43) = {fbxmin, fdymax, fdzmax, rd};
Point(44) = {fbxmax, fdymax, fdzmax, rd};

Point(45) = {fdxmin, fdymax, fbzmax, rd};
Point(46) = {fdxmax, fdymax, fbzmax, rd};
Point(47) = {fbxmin, fdymax, fbzmax, rd};
Point(48) = {fbxmax, fdymax, fbzmax, rd};

Point(49) = {fdxmin, fdymax, fbzmin, rd};
Point(50) = {fdxmax, fdymax, fbzmin, rd};
Point(51) = {fbxmin, fdymax, fbzmin, rd};
Point(52) = {fbxmax, fdymax, fbzmin, rd};

Point(53) = {fdxmin, fdymax, fdzmin, rd};
Point(54) = {fdxmax, fdymax, fdzmin, rd};
Point(55) = {fbxmin, fdymax, fdzmin, rd};
Point(56) = {fbxmax, fdymax, fdzmin, rd};


Line(192) = {54, 56};
Line(193) = {56, 52};
Line(194) = {52, 50};
Line(195) = {50, 54};
Line(196) = {50, 46};
Line(197) = {46, 48};
Line(198) = {48, 52};
Line(199) = {44, 48};
Line(200) = {42, 44};
Line(201) = {42, 46};
Line(202) = {35, 42};
Line(203) = {23, 46};
Line(204) = {19, 50};
Line(205) = {27, 54};
Line(206) = {44, 37};
Line(207) = {48, 15};
Line(208) = {52, 11};
Line(209) = {56, 29};
Line(210) = {48, 47};
Line(211) = {47, 43};
Line(212) = {44, 43};
Line(213) = {43, 39};
Line(214) = {47, 14};
Line(215) = {51, 47};
Line(216) = {52, 51};
Line(217) = {10, 51};
Line(218) = {51, 55};
Line(219) = {56, 55};
Line(220) = {55, 31};
Line(221) = {43, 41};
Line(222) = {41, 34};
Line(223) = {41, 45};
Line(224) = {45, 47};
Line(225) = {45, 22};
Line(226) = {51, 49};
Line(227) = {45, 49};
Line(228) = {49, 18};
Line(229) = {55, 53};
Line(230) = {53, 49};
Line(231) = {53, 26};
Line Loop(232) = {200, 206, -144, 202};
Plane Surface(233) = {232};
Line Loop(234) = {202, 201, -203, 141};
Plane Surface(235) = {234};
Line Loop(236) = {197, 207, 61, 203};
Plane Surface(237) = {236};
Line Loop(238) = {206, -148, -207, -199};
Plane Surface(239) = {238};
Line Loop(240) = {200, 199, -197, -201};
Plane Surface(241) = {240};
Surface Loop(242) = {241, 233, 239, 237, 235, 164};
Volume(243) = {242};
Line Loop(244) = {212, 213, 151, -206};
Plane Surface(245) = {244};
Line Loop(246) = {212, -211, -210, -199};
Plane Surface(247) = {246};
Line Loop(248) = {213, 152, -214, 211};
Plane Surface(249) = {248};
Line Loop(250) = {28, -214, -210, 207};
Plane Surface(251) = {250};
Surface Loop(252) = {247, 245, 249, 251, 239, 176};
Volume(253) = {252};
Line Loop(254) = {221, 222, 156, -213};
Plane Surface(255) = {254};
Line Loop(256) = {222, 157, -225, -223};
Plane Surface(257) = {256};
Line Loop(258) = {221, 223, 224, 211};
Plane Surface(259) = {258};
Line Loop(260) = {224, 214, 50, -225};
Plane Surface(261) = {260};
Surface Loop(262) = {259, 255, 257, 261, 185, 249};
Volume(263) = {262};


Line Loop(264) = {227, 228, -53, -225};
Plane Surface(265) = {264};
Line Loop(266) = {224, -215, 226, -227};
Plane Surface(267) = {266};
Line Loop(268) = {228, -57, 217, 226};
Plane Surface(269) = {268};
Line Loop(270) = {215, 214, 30, 217};
Plane Surface(271) = {270};
Surface Loop(272) = {267, 271, 269, 265, 67, 261};
Volume(273) = {272};
Line Loop(274) = {230, 228, 108, -231};
Plane Surface(275) = {274};
Line Loop(276) = {230, -226, 218, 229};
Plane Surface(277) = {276};
Line Loop(278) = {231, -104, -220, 229};
Plane Surface(279) = {278};
Line Loop(280) = {220, -103, 217, 218};
Plane Surface(281) = {280};
Surface Loop(282) = {277, 275, 279, 281, 269, 136};
Volume(283) = {282};
Line Loop(284) = {220, -99, -209, 219};
Plane Surface(285) = {284};
Line Loop(286) = {219, -218, -216, -193};
Plane Surface(287) = {286};
Line Loop(288) = {217, -216, 208, -32};
Plane Surface(289) = {288};
Line Loop(290) = {193, 208, -98, -209};
Plane Surface(291) = {290};
Surface Loop(292) = {287, 285, 291, 289, 281, 126};
Volume(293) = {292};
Line Loop(294) = {192, 209, -94, 205};
Plane Surface(295) = {294};
Plane Surface(296) = {294};
Line Loop(297) = {193, 194, 195, 192};
Plane Surface(298) = {297};
Line Loop(299) = {205, -195, -204, -93};
Plane Surface(300) = {299};
Line Loop(301) = {194, -204, 58, -208};
Plane Surface(302) = {301};
Surface Loop(303) = {298, 302, 300, 295, 118, 291};
Volume(304) = {303};
Line Loop(305) = {204, 196, -203, -62};
Plane Surface(306) = {305};
Line Loop(307) = {194, 196, 197, 198};
Plane Surface(308) = {307};
Line Loop(309) = {198, 208, 34, -207};
Plane Surface(310) = {309};
Surface Loop(311) = {308, 306, 310, 80, 237, 302};
Volume(312) = {311};
Line Loop(313) = {215, -210, 198, 216};
Plane Surface(314) = {313};
Surface Loop(315) = {314, 289, 310, 251, 271, 43};
Volume(316) = {315};


Coherence;
Physical Surface("Inflow", 2) = {275, 138, 71, 265, 189, 257, 134, 122, 112, 279, 295, 285, 162, 233, 245, 255, 183, 172, 298, 308, 241, 287, 277, 267, 314, 259, 247};
Physical Surface("Ourflow", 3) = {300, 110, 306, 86, 235, 160};
Physical Surface("Wall", 1) = {187, 75, 132, 174, 166, 84, 114, 37, 124, 20, 18, 16, 14, 22};
Physical Volume("Domain", 9) = {283, 293, 304, 312, 243, 120, 90, 170, 253, 316, 273, 263, 130, 49, 181, 191, 78, 140};
Coherence;
