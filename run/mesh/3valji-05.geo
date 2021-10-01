/*------------------------------------------------------------------

         File Name: pipe.geo

         Purpose:

         Creation Date: 05-04-2019

         Last Modified: sre 22 maj 2019 10:20:59 CEST

         Created By: teiresias

-----------------------------------------------------------------*/

L = 1; //        Length of the pipe
R = 0.5;  //      Diameter of a large pipe (no need to change)
lsize = 0.5; // element size

sqRR = R/Sqrt(2);

Point(0) = { 0, 0, 0, lsize }; 
Point(9) = { sqRR, sqRR, 0, lsize }; 
Point(10) = { sqRR, -sqRR, 0, lsize }; 
Point(11) = { -sqRR, -sqRR, 0, lsize }; 
Point(12) = { -sqRR, sqRR, 0, lsize }; 
//+
Circle(1) = {12, 0, 9};
//+
Circle(2) = {9, 0, 10};
//+
Circle(3) = {10, 0, 11};
//+
Circle(4) = {11, 0, 12};
//+
Line Loop(1) = {3, 4, 1, 2};
//+
Surface(1) = {1};
//+
Translate {0, 0, 1} {
  Duplicata { Point{12}; Point{9}; Point{0}; Point{11}; Point{10}; Line{1}; Line{4}; Line{2}; Line{3}; Surface{1}; }
}
//+
Rotate {{1, 0, 0}, {0, 0, 1}, Pi} {
  Duplicata { Point{12}; Point{9}; Point{0}; Point{11}; Point{10}; Line{1}; Line{4}; Line{2}; Line{3}; Surface{1}; }
}
Rotate {{1, 0, 0}, {0, 0, 1.5}, Pi} {
  Duplicata { Point{12}; Point{9}; Point{0}; Point{11}; Point{10}; Line{1}; Line{4}; Line{2}; Line{3}; Surface{1}; }
}
//+
Line(19) = {26, 21};
//+
Line(20) = {21, 13};
//+
Line(21) = {13, 12};
//+
Line(22) = {9, 14};
//+
Line(23) = {14, 22};
//+
Line(24) = {22, 27};
//+
Line(25) = {23, 18};
//+
Line(26) = {18, 16};
//+
Line(27) = {16, 11};
//+
Line(28) = {10, 17};
//+
Line(29) = {17, 19};
//+
Line(30) = {19, 24};
//+
Line Loop(2) = {22, 7, -28, -2};
//+
Surface(20) = {2};
//+
Line Loop(3) = {21, 1, 22, -5};
//+
Surface(21) = {-3};
//+
Line Loop(4) = {27, 4, -21, -6};
//+
Surface(22) = {-4};
//+
Line Loop(5) = {27, -3, 28, 8};
//+
Surface(23) = {5};
//+
Line Loop(6) = {6, -20, 11, 26};
//+
Surface(24) = {-6};
//+
Line Loop(7) = {20, 5, 23, 13};
//+
Surface(25) = {-7};
//+
Line Loop(8) = {26, -8, 29, -10};
//+
Surface(26) = {8};
//+
Line Loop(9) = {12, -23, 7, 29};
//+
Surface(27) = {-9};
//+
Line Loop(10) = {25, -11, -19, 16};
//+
Surface(28) = {-10};
//+
Line Loop(11) = {19, -13, 24, 18};
//+
Surface(29) = {-11};
//+
Line Loop(12) = {17, -24, -12, 30};
//+
Surface(30) = {-12};
//+
Line Loop(13) = {25, 10, 30, -15};
//+
Surface(31) = {13};
//+
Physical Surface("in1") = {1};
//+
Physical Surface("in3") = {19};
//+
Physical Surface("wall1") = {21, 22, 20, 23};
//+
Physical Surface("wall3") = {28, 29, 30, 31};
//+
Physical Surface("wall2") = {25, 24, 26, 27};
//+
Physical Surface("contact12") = {9};
//+
Physical Surface("contact23") = {14};
