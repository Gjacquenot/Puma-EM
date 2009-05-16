// some lengths

a_1 = 0.238; // first aperture
b_1 = 0.138; // first aperture
h_1 = 0.0; // depth of 'first aperture'

a_2 = 0.087; // outer contour of 'second aperture'
a_2b = 0.085; // inner contour of 'second aperture'
b_2 = 0.064; // inner contour of 'second aperture'
h_2 = -0.1555; // depth of 'second aperture'

a_3 = 0.025; // bottom aperture
b_3 = 0.016; // bottom aperture
h_3 = -0.189; // depth of 'bottom aperture'

w_ridge = 0.0075; // width of ridge

lc_coarse = 0.01;
lc_detail = 0.005;

// aperture

Point(1) = {-a_1/2.0, -b_1/2.0, h_1, lc_coarse};
Point(2) = {a_1/2.0, -b_1/2.0, h_1, lc_coarse};
Point(3) = {a_1/2.0, b_1/2.0, h_1, lc_coarse};
Point(4) = {-a_1/2.0, b_1/2.0, h_1, lc_coarse};


// outer contour of 'second aperture'

Point(5) = {-a_2/2.0, -b_2/2.0, h_2, lc_detail};
Point(6) = {a_2/2.0, -b_2/2.0, h_2, lc_detail};
Point(7) = {a_2/2.0, b_2/2.0, h_2, lc_detail};
Point(8) = {-a_2/2.0, b_2/2.0, h_2, lc_detail};


// inner contour of 'second aperture'

Point(9) = {-a_2b/2.0, -b_2/2.0, h_2, lc_detail};
Point(10) = {a_2b/2.0, -b_2/2.0, h_2, lc_detail};
Point(11) = {a_2b/2.0, b_2/2.0, h_2, lc_detail};
Point(12) = {-a_2b/2.0, b_2/2.0, h_2, lc_detail};


// contour of 'bottom aperture'

Point(13) = {-a_3/2.0, -b_3/2.0, h_3, lc_detail};
Point(14) = {a_3/2.0, -b_3/2.0, h_3, lc_detail};
Point(15) = {a_3/2.0, b_3/2.0, h_3, lc_detail};
Point(16) = {-a_3/2.0, b_3/2.0, h_3, lc_detail};


// 2 opposite slopes of 1.6 cm width from inner contour (shortest sides) of 'second aperture' to 'bottom aperture'

Point(17) = {-a_2b/2.0,-b_3/2.0, h_2,lc_detail};
Point(18) = {-a_2b/2.0,b_3/2.0, h_2,lc_detail};
Point(19) = {a_2b/2.0,-b_3/2.0, h_2,lc_detail};
Point(20) = {a_2b/2.0,b_3/2.0, h_2,lc_detail};
Point(21) = {-a_3/2.0,-b_3/2.0,-0.183,lc_detail};
Point(22) = {a_3/2.0,-b_3/2.0,-0.183,lc_detail};
Point(23) = {a_3/2.0,b_3/2.0,-0.183,lc_detail};
Point(24) = {-a_3/2.0,b_3/2.0,-0.183,lc_detail};


// 4 opposite slopes of 8.5 cm width from inner contour (longest sides) of 'second aperture' to 'bottom aperture'

Point(25) = {-a_2b/2.0, -b_3/2.0, -0.184, lc_detail};
Point(26) = {-a_2b/2.0, b_3/2.0, -0.184, lc_detail};
Point(27) = {a_2b/2.0, -b_3/2.0, -0.184, lc_detail};
Point(28) = {a_2b/2.0, b_3/2.0, -0.184, lc_detail};
Point(29) = {-a_3/2.0, -b_3/2.0, -0.184, lc_detail};
Point(30) = {a_3/2.0, -b_3/2.0, -0.184, lc_detail};
Point(31) = {a_3/2.0, b_3/2.0, -0.184, lc_detail};
Point(32) = {-a_3/2.0, b_3/2.0, -0.184, lc_detail};

// Ridge
// bottom points
Point(33) = {-w_ridge/2.0, -0.00075, -0.18, lc_detail};
Point(34) = {w_ridge/2.0, -0.00075, -0.18, lc_detail};
Point(35) = {w_ridge/2.0, 0.00075, -0.18, lc_detail};
Point(36) = {-w_ridge/2.0, 0.00075, -0.18, lc_detail};
Point(37) = {-w_ridge/2.0, -b_3/2.0, -0.18, lc_detail};
Point(38) = {w_ridge/2.0, -b_3/2.0, -0.18, lc_detail};
Point(39) = {w_ridge/2.0, b_3/2.0, -0.18, lc_detail};
Point(40) = {-w_ridge/2.0, b_3/2.0, -0.18, lc_detail};

// intersection with 'second aperture'
Point(41) = {-w_ridge/2.0, -b_2/2.0, h_2, lc_detail};
Point(42) = {w_ridge/2.0, -b_2/2.0, h_2, lc_detail};
Point(43) = {w_ridge/2.0, b_2/2.0, h_2, lc_detail};
Point(44) = {-w_ridge/2.0, b_2/2.0, h_2, lc_detail};

// intersection with 'first aperture'
Point(45) = {-w_ridge/2.0, -b_1/2.0, h_1, lc_coarse};
Point(46) = {w_ridge/2.0, -b_1/2.0, h_1, lc_coarse};
Point(47) = {w_ridge/2.0, b_1/2.0, h_1, lc_coarse};
Point(48) = {-w_ridge/2.0, b_1/2.0, h_1, lc_coarse};

// intersection with the 4 opposite slopes
Point(49) = {-w_ridge/2.0, -b_3/2.0, -0.184, lc_detail};
Point(50) = {w_ridge/2.0, -b_3/2.0, -0.184, lc_detail};
Point(51) = {w_ridge/2.0, b_3/2.0, -0.184, lc_detail};
Point(52) = {-w_ridge/2.0, b_3/2.0, -0.184, lc_detail};

// intersection with the 'bottom aperture'
Point(53) = {-w_ridge/2.0, -b_3/2.0, h_3, lc_detail};
Point(54) = {w_ridge/2.0, -b_3/2.0, h_3, lc_detail};
Point(55) = {w_ridge/2.0, b_3/2.0, h_3, lc_detail};
Point(56) = {-w_ridge/2.0, b_3/2.0, h_3, lc_detail};

// ridge profile (resulting from painful measurements) 
/*  Columns 1 through 6 (all in cm)
 *
 * -18.00000000000000 -10.10000000000000  -8.60000000000000  -7.40000000000000  -6.65000000000000  -6.10000000000000
 *   0.00075000000000   0.50000000000000   0.75000000000000   1.00000000000000   1.25000000000000   1.45000000000000
 *
 *  Columns 7 through 10
 *
 *  -4.90000000000000  -3.30000000000000  -2.10000000000000                  0
 *   2.00000000000000   3.00000000000000   4.20000000000000   6.90000000000000
 */


Point(57) = {-w_ridge/2.0, -0.0008, -0.17, lc_detail};
Point(58) = {w_ridge/2.0, -0.0008, -0.17, lc_detail};
Point(59) = {w_ridge/2.0, 0.0008, -0.17, lc_detail};
Point(60) = {-w_ridge/2.0, 0.0008, -0.17, lc_detail};

Point(61) = {-w_ridge/2.0, -0.005, -0.101, lc_detail};
Point(62) = {w_ridge/2.0, -0.005, -0.101, lc_detail};
Point(63) = {w_ridge/2.0, 0.005, -0.101, lc_detail};
Point(64) = {-w_ridge/2.0, 0.005, -0.101, lc_detail};

Point(65) = {-w_ridge/2.0, -0.0075, -0.086, lc_detail};
Point(66) = {w_ridge/2.0, -0.0075, -0.086, lc_detail};
Point(67) = {w_ridge/2.0, 0.0075, -0.086, lc_detail};
Point(68) = {-w_ridge/2.0, 0.0075, -0.086, lc_detail};

Point(69) = {-w_ridge/2.0, -0.01, -0.074, lc_detail};
Point(70) = {w_ridge/2.0, -0.01, -0.074, lc_detail};
Point(71) = {w_ridge/2.0, 0.01, -0.074, lc_detail};
Point(72) = {-w_ridge/2.0, 0.01, -0.074, lc_detail};

Point(73) = {-w_ridge/2.0, -0.0125, -0.0665, lc_detail};
Point(74) = {w_ridge/2.0, -0.0125, -0.0665, lc_detail};
Point(75) = {w_ridge/2.0, 0.0125, -0.0665, lc_detail};
Point(76) = {-w_ridge/2.0, 0.0125, -0.0665, lc_detail};

Point(77) = {-w_ridge/2.0, -0.0145, -0.061, lc_detail};
Point(78) = {w_ridge/2.0, -0.0145, -0.061, lc_detail};
Point(79) = {w_ridge/2.0, 0.0145, -0.061, lc_detail};
Point(80) = {-w_ridge/2.0, 0.0145, -0.061, lc_detail};

Point(81) = {-w_ridge/2.0, -0.02, -0.049, lc_detail};
Point(82) = {w_ridge/2.0, -0.02, -0.049, lc_detail};
Point(83) = {w_ridge/2.0, 0.02, -0.049, lc_detail};
Point(84) = {-w_ridge/2.0, 0.02, -0.049, lc_detail};

Point(85) = {-w_ridge/2.0, -0.03, -0.033, lc_detail};
Point(86) = {w_ridge/2.0, -0.03, -0.033, lc_detail};
Point(87) = {w_ridge/2.0, 0.03, -0.033, lc_detail};
Point(88) = {-w_ridge/2.0, 0.03, -0.033, lc_detail};

Point(89) = {-w_ridge/2.0, -0.042, -0.021, lc_detail};
Point(90) = {w_ridge/2.0, -0.042, -0.021, lc_detail};
Point(91) = {w_ridge/2.0, 0.042, -0.021, lc_detail};
Point(92) = {-w_ridge/2.0, 0.042, -0.021, lc_detail};

// lines

// not classified



Line(1) = {1,45};
Line(2) = {45,46};
Line(3) = {46,2};
Line(4) = {2,3};
Line(5) = {3,47};
Line(6) = {47,48};
Line(7) = {48,4};
Line(8) = {4,1};
Line(9) = {9,41};
Line(10) = {41,42};
Line(11) = {42,10};
Line(12) = {6,7};
Line(13) = {11,43};
Line(14) = {43,44};
Line(15) = {44,12};
Line(16) = {8,5};
Line(17) = {12,18};
Line(18) = {18,17};
Line(19) = {17,9};
Line(20) = {11,20};
Line(21) = {20,19};
Line(22) = {19,10};
Line(23) = {10,6};
Line(24) = {7,11};
Line(25) = {12,8};
Line(26) = {5,9};
Line(27) = {4,8};
Line(28) = {3,7};
Line(29) = {6,2};
Line(30) = {1,5};
Line(31) = {45,89};
Line(32) = {89,85};
Line(33) = {85,81};
Line(34) = {81,77};
Line(35) = {77,73};
Line(36) = {73,69};
Line(37) = {69,65};
Line(38) = {65,61};
Line(39) = {62,66};
Line(40) = {66,70};
Line(41) = {70,74};
Line(42) = {74,78};
Line(43) = {78,82};
Line(44) = {82,86};
Line(45) = {86,90};
Line(46) = {90,46};
Line(47) = {48,92};
Line(48) = {92,88};
Line(49) = {88,84};
Line(50) = {84,80};
Line(51) = {80,76};
Line(52) = {76,72};
Line(53) = {72,68};
Line(54) = {68,64};
Line(55) = {47,91};
Line(56) = {91,87};
Line(57) = {87,83};
Line(58) = {83,79};
Line(59) = {79,75};
Line(60) = {75,71};
Line(61) = {71,67};
Line(62) = {67,63};
Line(63) = {43,47};
Line(64) = {48,44};
Line(65) = {46,42};
Line(66) = {45,41};
Line(67) = {18,26};
Line(68) = {26,12};
Line(69) = {25,17};
Line(70) = {25,9};
Line(71) = {20,28};
Line(72) = {28,11};
Line(73) = {19,27};
Line(74) = {27,10};
Line(75) = {26,32};
Line(76) = {32,52};
Line(77) = {52,51};
Line(78) = {51,31};
Line(79) = {31,28};
Line(80) = {25,29};
Line(81) = {29,49};
Line(82) = {49,50};
Line(83) = {50,30};
Line(84) = {30,27};
Line(85) = {13,53};
Line(86) = {53,54};
Line(87) = {54,14};
Line(88) = {14,15};
Line(89) = {15,55};
Line(90) = {55,56};
Line(91) = {56,16};
Line(92) = {16,13};
Line(93) = {52,44};
Line(94) = {43,51};
Line(95) = {49,41};
Line(96) = {42,50};
Line(97) = {32,24};
Line(98) = {24,18};
Line(99) = {17,21};
Line(100) = {21,29};
Line(101) = {30,22};
Line(102) = {22,19};
Line(103) = {31,23};
Line(104) = {23,20};
Line(105) = {52,40};
Line(106) = {40,36};
Line(107) = {33,37};
Line(108) = {37,49};
Line(109) = {50,38};
Line(110) = {38,34};
Line(111) = {35,39};
Line(112) = {39,51};
Line(113) = {36,60};
Line(114) = {33,57};
Line(115) = {35,59};
Line(116) = {34,58};
Line(117) = {57,61};
Line(118) = {58,62};
Line(119) = {64,60};
Line(120) = {59,63};
Line(121) = {16,32};
Line(122) = {13,29};
Line(123) = {14,30};
Line(124) = {15,31};
Line(137) = {92,91};
Line(138) = {88,87};
Line(139) = {84,83};
Line(140) = {80,79};
Line(141) = {76,75};
Line(142) = {72,71};
Line(143) = {68,67};
Line(144) = {64,63};
Line(145) = {60,59};
Line(146) = {58,57};
Line(147) = {36,35};
Line(148) = {34,33};
Line(149) = {62,61};
Line(150) = {66,65};
Line(151) = {70,69};
Line(152) = {74,73};
Line(153) = {78,77};
Line(154) = {82,81};
Line(155) = {86,85};
Line(156) = {90,89};
Line(181) = {39,40};
Line(182) = {38,37};
Line(233) = {22,23};
Line(234) = {21,24};


Line Loop(125) = {11,23,29,-3,65};
Plane Surface(126) = {125};
Line Loop(127) = {10,-65,-2,66};
Plane Surface(128) = {127};
Line Loop(129) = {66,-9,-26,-30,1};
Plane Surface(130) = {129};
Line Loop(131) = {28,24,13,63,-5};
Plane Surface(132) = {131};
Line Loop(133) = {14,-64,-6,-63};
Plane Surface(134) = {133};
Line Loop(135) = {7,27,-25,-15,-64};
Plane Surface(136) = {135};
Line Loop(157) = {64,-93,105,106,113,-119,-54,-53,-52,-51,-50,-49,-48,-47};
Plane Surface(158) = {157};
Line Loop(159) = {63,55,56,57,58,59,60,61,62,-120,-115,111,112,-94};
Plane Surface(160) = {159};
Line Loop(161) = {55,-137,-47,-6};
Plane Surface(162) = {161};
Line Loop(163) = {48,138,-56,-137};
Plane Surface(164) = {163};
Line Loop(165) = {57,-139,-49,138};
Plane Surface(166) = {165};
Line Loop(167) = {58,-140,-50,139};
Plane Surface(168) = {167};
Line Loop(169) = {59,-141,-51,140};
Plane Surface(170) = {169};
Line Loop(171) = {142,-60,-141,52};
Plane Surface(172) = {171};
Line Loop(173) = {143,-61,-142,53};
Plane Surface(174) = {173};
Line Loop(175) = {144,-62,-143,54};
Plane Surface(176) = {175};
Line Loop(177) = {145,120,-144,119};
Plane Surface(178) = {177};
Line Loop(179) = {115,-145,-113,147};
Plane Surface(180) = {179};
Line Loop(183) = {147,111,181,106};
Plane Surface(184) = {183};
Line Loop(185) = {112,-77,105,-181};
Plane Surface(186) = {185};
Line Loop(187) = {77,-94,14,-93};
Plane Surface(188) = {187};
Line Loop(189) = {2,-46,156,-31};
Plane Surface(190) = {189};
Line Loop(191) = {45,156,32,-155};
Plane Surface(192) = {191};
Line Loop(193) = {44,155,33,-154};
Plane Surface(194) = {193};
Line Loop(195) = {43,154,34,-153};
Plane Surface(196) = {195};
Line Loop(197) = {42,153,35,-152};
Plane Surface(198) = {197};
Line Loop(199) = {41,152,36,-151};
Plane Surface(200) = {199};
Line Loop(201) = {40,151,37,-150};
Plane Surface(202) = {201};
Line Loop(203) = {39,150,38,-149};
Plane Surface(204) = {203};
Line Loop(205) = {118,149,-117,-146};
Plane Surface(206) = {205};
Line Loop(207) = {148,114,-146,-116};
Plane Surface(208) = {207};
Line Loop(209) = {110,148,107,-182};
Plane Surface(210) = {209};
Line Loop(211) = {82,109,182,108};
Plane Surface(212) = {211};
Line Loop(213) = {96,-82,95,10};
Plane Surface(214) = {213};
Line Loop(215) = {65,96,109,110,116,118,39,40,41,42,43,44,45,46};
Plane Surface(216) = {215};
Line Loop(217) = {95,-66,31,32,33,34,35,36,37,38,-117,-114,107,108};
Plane Surface(218) = {217};
Line Loop(219) = {27,16,-30,-8};
Plane Surface(220) = {219};
Line Loop(221) = {12,-28,-4,-29};
Plane Surface(222) = {221};
Line Loop(223) = {86,87,88,89,90,91,92,85};
Plane Surface(224) = {223};
Line Loop(225) = {72,20,71};
Plane Surface(226) = {225};
Line Loop(227) = {73,74,-22};
Plane Surface(228) = {227};
Line Loop(229) = {68,17,67};
Plane Surface(230) = {229};
Line Loop(231) = {69,19,-70};
Plane Surface(232) = {231};
Line Loop(235) = {234,98,18,99};
Plane Surface(236) = {235};
Line Loop(237) = {102,-21,-104,-233};
Plane Surface(238) = {237};
Line Loop(239) = {233,-103,-124,-88,123,101};
Plane Surface(240) = {239};
Line Loop(241) = {121,97,-234,100,-122,-92};
Plane Surface(242) = {241};
Line Loop(243) = {123,-83,-82,-81,-122,85,86,87};
Plane Surface(244) = {243};
Line Loop(245) = {76,77,78,-124,89,90,91,121};
Plane Surface(246) = {245};
Line Loop(247) = {99,100,-80,69};
Plane Surface(248) = {247};
Line Loop(249) = {98,67,75,97};
Plane Surface(250) = {249};
Line Loop(251) = {102,73,-84,101};
Plane Surface(252) = {251};
Line Loop(253) = {19,-26,-16,-25,17,18};
Plane Surface(254) = {253};
Line Loop(255) = {12,24,20,21,22,23};
Plane Surface(256) = {255};
Line Loop(257) = {80,81,95,-9,-70};
Plane Surface(258) = {257};
Line Loop(259) = {11,-74,-84,-83,-96};
Plane Surface(260) = {259};
Line Loop(261) = {93,15,-68,75,76};
Plane Surface(262) = {261};
Line Loop(263) = {78,79,72,13,94};
Plane Surface(264) = {263};
Line Loop(265) = {104,71,-79,103};
Plane Surface(266) = {265};
