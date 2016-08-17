clear all; close all; clc;
%% Implementation
% Knot vector
Xi = [ 0 0 0 1 2 3 4 4 5 5 5];

% Coordinates for B-spline/NURBS curve
x = [0 1  2  2 4 5 2 1]';
y = [0 -1 -1 1 1 3 5 2]' * .7;
z = [0  0  0 0 0 0 0 0]';

% Weights (if all equal -> B-spline, else NURBS);
w = [1 1 2 1 1 1 3 1]';

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

%% Generate basis
res = 40; % Compute 40 samples, uniformly distributed (and knot points are auto-included)
[ R , U] = nrbasis_num( Curve.KV, Curve.w, res );

%% Calculate NURBS-spline curve
C = R' * Curve.CP;

%% Plot basis
figure(1)
plotNurbsBasis(R,U);

%% PLot curve
figure()
plotNurbsCurve2D(C(:,1),C(:,2),Curve.KV,U,x,y);