clear all; close all; clc;
%% Read data
% Load data from a (.m) file
tic
% data_solid
% data_solid2
% data_solid3
data_solid4
% data_solid5
% data_solid6
% data_solid7
% data_solid8
% data_solid9

fprintf('Time to load data: %2.2f seconds\n',toc); tic;
%% Perform knot refinement
% Refine once in Xi
[ X_ ] = getSubDivKVValues( Xi, 1 );
[B,Xi] = RefineKnotSolidXi( B,p,Xi,X_ );

% Refine twice in Eta and Zeta
[ E_ ] = getSubDivKVValues( Eta, 2 );
[ B,Eta ] = RefineKnotSolidEta( B,q,Eta,E_ );
[ Z_ ] = getSubDivKVValues( Zeta, 2 );
[ B,Zeta ] = RefineKnotSolidZeta( B,r,Zeta,Z_ );

%% Compute basis functions and mapping
res_x = 20;
res_y = 20;
res_z = 20;

% Number of non-zero elements:
nxi=0;
neta=0;
nzeta=0;
for i = 1 : length(Xi)-1
    if (Xi(i)) ~= Xi(i+1)
        nxi = nxi+1;
    end
end
for i = 1 : length(Eta)-1
    if (Eta(i)) ~= Eta(i+1)
        neta = neta+1;
    end
end
for i = 1 : length(Zeta)-1
    if (Zeta(i)) ~= Zeta(i+1)
        nzeta = nzeta+1;
    end
end
fprintf('Number of nnz elements: %i \n',neta*nxi*nzeta)

tic
% Nurbs surface basis from given knots and weigths(in Control point data).
[R, U, V, W] = nrbasis_solid_num(Xi,Eta,Zeta,B,res_x,res_y,res_z);
fprintf('Time to generate basis: %2.2f seconds\n',toc); tic;

% Nurbs surface according to (2.29) in Cotrell, Hughes & Bazilevs
Cx = zeros(size(R{1}));
Cy = zeros(size(R{1}));
Cz = zeros(size(R{1}));
for i = 1:size(R,1)
    for j = 1 : size(R,2)
        for k = 1 : size(R,3)
        Cx = Cx + R{i,j,k}*B{i,j,k}(1);
        Cy = Cy + R{i,j,k}*B{i,j,k}(2);
        Cz = Cz + R{i,j,k}*B{i,j,k}(3);
        end
    end
end
fprintf('Time to generate surface: %2.2f seconds\n',toc); tic;

%% Plotting the Nurbs solid
figure(1)
plotNurbsSolid( Cx,Cy,Cz, B )
fprintf('Time to plot surface: %2.2f seconds\n',toc); tic;
% Plotting the knot lines
figure(2)
plotNurbsSolidElement( Cx, Cy, Cz, Xi, Eta, Zeta, U, V, W )
fprintf('Time to plot elements: %2.2f seconds\n',toc); tic;