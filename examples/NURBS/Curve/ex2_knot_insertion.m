clear all; close all; clc;
%%
resolution = 100; % Number of places to evaluate the curve at

% Knot vector
Xi = [ 0 0 0 1 2 3 4 4 5 5 5];

% Control points:
x = [0 1  2  2 4 5 2 1]';
y = [0 -1 -1 1 1 3 5 2]';
w = [1  1  1 2 3 1 5 1]';
% 
% Number of knots
k = length(Xi);

% Number of control-points
n = length(x);

% Order of basis
p = k-n-1;

% Calculate NURBS basis
[ R , Xi_store] = nrbasis_num( Xi, w, resolution );

% Plot knot-vector vs NURBS basis functions
subplot(2,2,1)
plot(Xi_store,R)
grid on
title('NURBS Basis functions')
xlabel('$\xi$','interpreter','latex')

% Plot Curve (NURBS)
% Calculate the NURBS curve.
C2 = [x,y]' * R;

% Plot:
subplot(2,2,2)
plot(C2(1,:),C2(2,:))
hold on
title('NURBS curve')
xlabel('x')
ylabel('y')

% Plot the control-points
plot(x,y,'ro--')
grid on

% Plot the knots:
knots = unique(Xi);
knot_pos = zeros(length(knots),1);
for i = 1 : length(knots)
    knot_pos(i) = find(Xi_store == knots(i)); % Find the knot position
end
plot(C2(1,knot_pos),C2(2,knot_pos),'bs')
legend('Curve','Control Points','Knots')
axis equal

%% Do knot insertion

% Knots to insert:
knot_insert = [ 0.315 ]; 

Pw = [x.*w y.*w w];            % We do the refinement by treating the 2D NURBS
r = length(knot_insert)-1;     % as a 3D B-Spline (see part 2.2.1 in IGA book)
                               % and use the B-Spline knot-insertion function:
[Xi_bar,Qw] = RefineKnotVectCurve(length(x)-1,p,Xi,Pw,knot_insert,r);
qw = Qw(:,3);
Qw(:,1) = Qw(:,1)./qw;
Qw(:,2) = Qw(:,2)./qw;


%% Calculate new basis:
Xi = Xi_bar; % Replace old knot-vector
x = Qw(:,1); y = Qw(:,2); w = Qw(:,3); % Replace control points
% Number of knots
k = length(Xi_bar);

% Number of control-points
n = length(x);

% Order of basis
p = k-n-1;

% Calculate NURBS basis
[ R , Xi_store] = nrbasis_num( Xi, w, resolution );

% Plot knot-vector vs NURBS basis functions
subplot(2,2,3)
plot(Xi_store,R)
grid on
title('NURBS Basis functions after knot insertion')
xlabel('$\xi$','interpreter','latex')

% Plot Curve (NURBS)
% Calculate the NURBS curve.
C2 = [x,y]' * R;

% Plot:
subplot(2,2,4)
plot(C2(1,:),C2(2,:))
hold on
title('NURBS curve after knot insertion')
xlabel('x')
ylabel('y')

% Plot the control-points
plot(x,y,'ro--')
grid on

% Plot the knots:
knots = unique(Xi);
knot_pos = zeros(length(knots),1);
for i = 1 : length(knots)
    knot_pos(i) = find(Xi_store == knots(i)); % Find the knot position
end
plot(C2(1,knot_pos),C2(2,knot_pos),'bs')
legend('Curve','Control Points','Knots')
axis equal