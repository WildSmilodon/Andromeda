/*------------------------------------------------------------------

         File Name: surf_hexa-kocka.geo

         Purpose:

         Creation Date: 22-05-2019

         Last Modified: sre 22 maj 2019 12:04:49 CEST

         Created By: teiresias

-----------------------------------------------------------------*/
//+
SetFactory("OpenCASCADE");
Box(1) = {0, -0, 0, 1, 1, 1};
//+
Physical Surface("spredaj") = {6};
//+
Physical Surface("zadaj") = {5};
//+
Physical Surface("zgoraj") = {4};
//+
Physical Surface("spodaj") = {3};
//+
Physical Surface("levo") = {1};
//+
Physical Surface("desno") = {2};
//+
Point(9) = {0.4, 1.2, 0, 1.0};
//+
Point(10) = {0.4, 1.2, 1, 1.0};
//+
Point(11) = {0.4, 0.5, 0, 1.0};
//+
Point(12) = {0.4, 0.5, 1, 1.0};
