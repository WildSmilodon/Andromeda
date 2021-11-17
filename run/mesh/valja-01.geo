/*------------------------------------------------------------------

         File Name: pipe.geo

         Purpose:

         Creation Date: 05-04-2019

         Last Modified: sre 08 maj 2019 13:55:40 CEST

         Created By: teiresias

-----------------------------------------------------------------*/

L = 1; //        Length of the pipe
R = 0.5;  //      Diameter 

sqRR = R/Sqrt(2);
lsize = 0.1;

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
//+
Line(14) = {10, 17};
//+
Line(15) = {17, 19};
//+
Line(16) = {9, 14};
//+
Line(17) = {14, 22};
//+
Line(18) = {12, 13};
//+
Line(19) = {13, 21};
//+
Line(20) = {11, 16};
//+
Line(21) = {16, 18};
//+
Line Loop(2) = {16, 7, -14, -2};
//+
Surface(15) = {2};
//+
Line Loop(3) = {1, 16, -5, -18};
//+
Surface(16) = {-3};
//+
Line Loop(4) = {14, 8, -20, -3};
//+
Surface(17) = {4};
//+
Line Loop(5) = {6, -18, -4, 20};
//+
Surface(18) = {5};
//+
Line Loop(6) = {8, 21, 10, -15};
//+
Surface(19) = {-6};
//+
Line Loop(7) = {15, 12, -17, 7};
//+
Surface(20) = {-7};
//+
Line Loop(8) = {13, -19, 5, 17};
//+
Surface(21) = {-8};
//+
Line Loop(9) = {21, -11, -19, -6};
//+
Surface(22) = {9};
//+
//+
Surface Loop(1) = {16, 18, 1, 17, 15, 9};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {9, 21, 20, 14, 22, 19};
//+
Volume(2) = {2};
Physical Surface("in1") = {1};
//+
Physical Surface("in2") = {14};
//+
Physical Surface("wall1") = {18, 16, 17, 15};
//+
Physical Surface("wall2") = {19, 22, 21, 20};
//+
Physical Surface("contact") = {9};
