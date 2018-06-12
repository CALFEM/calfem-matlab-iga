clear all; close all; clc;
%% Implementation
% Knot vector
Xi=[0 0 0 0.25 0.25 0.5 0.5 0.75 0.75 1 1 1];
% Coordinates for B-spline/NURBS curve
R = 3;
x=[1 1 0 -1 -1 -1 0 1 1]*R;
x = x';
y=[0 1 1 1 0 -1 -1 -1 0]*R;
y=y';
z=[0 0 0 0 0 0 0 0 0]';
% Weights (if all equal -> B-spline, else NURBS);
W = cos(pi/4);
w=[1 W 1 W 1 W 1 W 1]';




% Number of knots
k = length(Xi);
% Number of coordinates
n = length(x);
% Order of basis
p = k-n-1;

Curve.KV = Xi;
Curve.CP = [x,y,z];
Curve.w = w;
Curve.degree = p;
Curve.nCP = n;

%% Do plot
plotNurbsCurve2DSimple( Curve )
