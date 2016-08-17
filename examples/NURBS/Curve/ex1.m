clear all; close all; clc;
%%
resolution = 15; % Number of places to evaluate the curve at

% Knot vector
% Xi = [ 0 0 0 0 2 5 5 5 5];
Xi = [ 0 0 0 0 3 5 5 5 5];

% Control points:
x = [0 1  2    2 1]';
y = [0 -1 -1   5 2]';

% Weights:
w = [1  1  1   5 1]';

% Number of knots
k = length(Xi);

% Number of control-points
n = length(x);

% Order of basis
p = k-n-1;

%% Plot Basis
Xi_store =  min(Xi) : max(Xi)/resolution : max(Xi); % Evaluate basis functions at each xi in Xi_store
Xi_store = unique(sort([Xi_store Xi])); % <- Makes sure that the knots are included (In order to plot the knots, otherwise we could ignore this)
N = zeros(n,length(Xi_store)); % We have n basisfunctions, evaluated at each xi in Xi_store
for i = 1:length(Xi_store)
    xi = Xi_store(i); % Get current xi
    j=FindSpan(n,p,xi,Xi); % Find knot span
    N_loc=BasisFun(j,xi,p,Xi); % Get only the non-zero basis functions (there are p + 1 non-zero basis functions).
    N(j-p+1:j+1,i)=N_loc; % Store the non-zero basis functions
end

% Plot knot-vector vs B-Spline basis functions
subplot(2,2,1)
plot(Xi_store,N)
grid on
title('B-Spline Basis functions')
xlabel('$\xi$','interpreter','latex')

%% Plot Curve (B-spline)
% Calculate the curve. 
C = [x,y]' * N;

% Plot the curve:
subplot(2,2,2)
plot(C(1,:),C(2,:))
hold on
title('B-Spline curve')
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
plot(C(1,knot_pos),C(2,knot_pos),'bs')
legend('Curve','Control Points','Knots')
axis equal
%% Plot NURBS basis
% Calculate the NURBS basis
W = w'*N;  % Nurbs denominator (Sum of weights * basis functions)
R = zeros(size(N));
for i = 1 : size(R,2)
    R(:,i) = N(:,i).*w/W(i); % Calculate NURBS basis from B-Spline basis
end

% Plot knot-vector vs NURBS basis functions
subplot(2,2,3)
plot(Xi_store,R)
grid on
title('NURBS Basis functions')
xlabel('$\xi$','interpreter','latex')


%% Plot Curve (NURBS)
% Calculate the NURBS curve.
C2 = [x,y]' * R;

% Plot:
subplot(2,2,4)
plot(C2(1,:),C2(2,:))
hold on
title('NURBS curve')
xlabel('x')
ylabel('y')

% Plot the control-points
plot(x,y,'ro--')
grid on

% Plot the knots:
plot(C2(1,knot_pos),C2(2,knot_pos),'bs')
legend('Curve','Control Points','Knots')
axis equal