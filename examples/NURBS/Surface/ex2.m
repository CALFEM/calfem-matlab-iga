clear all; close all; clc;
% Surface example, consisting of "2 isogeometric elements"

resolution_xi = 20; % Number of places to evaluate the surface at in xi
resolution_eta = 20; % Number of places to evaluate the surface at in eta

% Knot vectors
Xi = [ 0 0 0  2 5 5 5];
Eta = [ 0 0 0 5 5 5];

% Control point net:
B=cell(4,3);
%        [x y z    w]
B{1,1} = [0 0 0    1];
B{2,1} = [1 0 1/2  1];
B{3,1} = [2 0 2/2  1];
B{4,1} = [3 0 3/2  1];
B{1,2} = [0 1 0    1];
B{2,2} = [1 1 1/2  1];
B{3,2} = [2 1 1/2  1];
B{4,2} = [3 1 1/2  5]; % w ~= 1
B{1,3} = [0 2 1/2  1];
B{2,3} = [1 2 2/2  1];
B{3,3} = [2 2 2/2  1];
B{4,3} = [3 2 2/2  1];



% Number of knots
k = length(Xi);
l = length(Eta);
% Number of control-points
n = size(B,1);
m = size(B,2);
% Order of basis
p = k-n-1;
q = l-m-1;

% Convert cells to matrices
for i = 1 : size(B,1)
    for j = 1 : size(B,2)
        B1(:,i,j) = B{i,j}(:); 
    end
end

Xi_store =  min(Xi) : max(Xi)/resolution_xi : max(Xi); % Evaluate basis functions at each xi in Xi_store
Xi_store = unique(sort([Xi_store Xi])); % <- Makes sure that the knots are included (In order to plot the knots, otherwise we could ignore this)
Eta_store =  min(Eta) : max(Eta)/resolution_xi : max(Eta); % Evaluate basis functions at each eta in Eta_store
Eta_store= unique(sort([Eta_store Eta])); % <- Makes sure that the knots are included (In order to plot the knots)

for i = 1:length(Xi_store)
    for j = 1:length(Eta_store)
        
        xi = Xi_store(i); % Get current xi
        eta = Eta_store(j); % Get current eta

        % Find xi span and get b-spline basis
        xiSpan = FindSpan(n,p,xi,Xi); % Find knot span
        N=BasisFun(xiSpan,xi,p,Xi);  % Get only the non-zero basis functions (there are p + 1 non-zero basis functions).
        
        % Find eta span and get b-spline basis
        etaSpan = FindSpan(m,q,eta,Eta); % Find knot span
        M=BasisFun(etaSpan,eta,q,Eta);  % Get only the non-zero basis functions (there are q + 1 non-zero basis functions).
        
        % Construct the bi-variate b-spline basis using the outer product
        NM=N'*M;
        
        % Get relevant control points for current non-zero basis funcitons.
        X=reshape(B1(1,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);
        Y=reshape(B1(2,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);
        Z=reshape(B1(3,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);
        w=reshape(B1(4,xiSpan-p+1:xiSpan+1,etaSpan-q+1:etaSpan+1),p+1,q+1);
        
        R = NM .* w; % Multiply Bspline basis by weights
        R = R / sum(sum(R)); % Divide by total W -> NURBS basis
        
        % Current point
        x = sum(sum(X.*R));
        y = sum(sum(Y.*R));
        z = sum(sum(Z.*R));
        
        % Save current surface point
        Sx(i,j) = x;
        Sy(i,j) = y;
        Sz(i,j) = z;
        
    end
end

% Plot the surface:
surf(Sx,Sy,Sz)
shading interp
colormap summer
hidden off
axis equal
title('NURBS surface')
xlabel('x')
ylabel('y')
zlabel('z')

% Plot the control points:
hold on
plot3(B1(1,:),B1(2,:),B1(3,:),'ro')

% To plot the knots we must find their positions in Xi_store, Eta_store
% Find knot positions
Xiu=unique(Xi);
xi_pos=[];
for i = 1 : length(Xiu)
    xi_pos = [xi_pos find(Xi_store == Xiu(i))];
end

% Find knot positions
Etau=unique(Eta);
eta_pos=[];
for i = 1 : length(Etau)
    eta_pos = [eta_pos find(Eta_store == Etau(i))];
end

% Plot the knots:
plot3(Sx(xi_pos,:)',Sy(xi_pos,:)',Sz(xi_pos,:)','b','LineWidth',2)
plot3(Sx(:,eta_pos),Sy(:,eta_pos),Sz(:,eta_pos),'b','LineWidth',2)

