//
//  gmsh mreza za kroglo v krogli
//

// Novi gmsh ne zasuka normal, uporabi ali
// gmsh -2 -format msh2 sphere-in-channel.geo
// ali
// ukaz ReverseMesh Surface{16}; to naredi



s=0.1; // Velikost elementa v krogli
Mesh.CharacteristicLengthFromCurvature = 0.1;
r=0.5; // radij krogle

// Velika krogla
ss=150.0; // Velikost elementa v krogli
rr=512.0; // radij krogle



// ********************************************************************* 2972-5936
// **************** MALA KROGLA *********************************************
// *********************************************************************

Point(1) = {0, 0, 0, s};
Point(2) = {-r, 0, 0, s};
Point(3) = {r, 0, 0, s};
Point(4) = {0, -r, 0, s};
Point(5) = {0, r, 0, s};
Point(6) = {0, 0, r, s};
Point(7) = {0, 0, -r, s};
Circle(1) = {3, 1, 5};
Circle(2) = {2, 1, 5};
Circle(3) = {2, 1, 4};
Circle(4) = {4, 1, 3};
Circle(5) = {6, 1, 3};
Circle(6) = {6, 1, 2};
Circle(7) = {2, 1, 7};
Circle(8) = {7, 1, 3};
Circle(9) = {5, 1, 6};
Circle(10) = {6, 1, 4};
Circle(11) = {4, 1, 7};
Circle(12) = {7, 1, 5};
Line Loop(13) = {9, 5, 1};
Surface(14) = {13};
Line Loop(15) = {2, 9, 6};
Surface(16) = {15};
Line Loop(17) = {7, 12, -2};
Surface(18) = {17};
Line Loop(19) = {12, -1, -8};
Surface(20) = {19};
Line Loop(21) = {6, 3, -10};
Surface(22) = {21};
Line Loop(23) = {10, 4, -5};
Surface(24) = {23};
Line Loop(25) = {11, -7, 3};
Surface(26) = {25};
Line Loop(27) = {8, -4, 11};
Surface(28) = {27};

// CASE : NORMALS OUT OF SPHERE
// Reverse normals
// ReverseMesh Surface{16}; // this was done using "-" in Physical Sruface in v2 mesh, in v4 it is not
// ReverseMesh Surface{18}; // this was automatic in v2 mesh, in v4 it is not
// ReverseMesh Surface{26}; // this was automatic in v2 mesh, in v4 it is not

// CASE : NORMALS INTO THE SPHERE
// Reverse normals
ReverseMesh Surface{22};
ReverseMesh Surface{20};
ReverseMesh Surface{14};
ReverseMesh Surface{24};
ReverseMesh Surface{28};

//Physical surfaces of domain
Physical Surface("particle") = {14, 16, 18, 20, 22, 24, 26, 28};




// *********************************************************************
// **************** MALA KROGLA *********************************************
// *********************************************************************

Point(11) = {0, 0, 0, ss};
Point(12) = {-rr, 0, 0, ss};
Point(13) = {rr, 0, 0, ss};
Point(14) = {0, -rr, 0, ss};
Point(15) = {0, rr, 0, ss};
Point(16) = {0, 0, rr, ss};
Point(17) = {0, 0, -rr, ss};
Circle(101) = {13, 11, 15};
Circle(102) = {12, 11, 15};
Circle(103) = {12, 11, 14};
Circle(104) = {14, 11, 13};
Circle(105) = {16, 11, 13};
Circle(106) = {16, 11, 12};
Circle(107) = {12, 11, 17};
Circle(108) = {17, 11, 13};
Circle(109) = {15, 11, 16};
Circle(110) = {16, 11, 14};
Circle(111) = {14, 11, 17};
Circle(112) = {17, 11, 15};
Line Loop(113) = {109, 105, 101};
Surface(114) = {113};
Line Loop(115) = {102, 109, 106};
Surface(116) = {115};
Line Loop(117) = {107, 112, -102};
Surface(118) = {117};
Line Loop(119) = {112, -101, -108};
Surface(120) = {119};
Line Loop(121) = {106, 103, -110};
Surface(122) = {121};
Line Loop(123) = {110, 104, -105};
Surface(124) = {123};
Line Loop(125) = {111, -107, 103};
Surface(126) = {125};
Line Loop(127) = {108, -104, 111};
Surface(128) = {127};

// CASE : NORMALS OUT OF SPHERE
// Reverse normals
ReverseMesh Surface{116}; // this was done using "-" in Physical Sruface in v2 mesh, in v4 it is not
ReverseMesh Surface{118}; // this was automatic in v2 mesh, in v4 it is not
ReverseMesh Surface{126}; // this was automatic in v2 mesh, in v4 it is not

// // CASE : NORMALS INTO THE SPHERE
// // Reverse normals
// ReverseMesh Surface{22};
// ReverseMesh Surface{20};
// ReverseMesh Surface{14};
// ReverseMesh Surface{24};
// ReverseMesh Surface{28};

//Physical surfaces of domain
Physical Surface("outside") = {114, 116, 118, 120, 122, 124, 126, 128};



