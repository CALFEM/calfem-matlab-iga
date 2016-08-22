clear all; close all; clc;

%% Load data
count=1;
run ./../NURBS/solid/data_solid6
% run ./../NURBS/solid/data_solid7
deg.p = p; deg.q=q; deg.r = r; clear p q r

% Perform 3 subdivisons in Xi, and 1 in  Eta & Zeta
[ X_ ] = getSubDivKVValues( Xi, 3 );
[B,Xi]=RefineKnotSolidXi( B,deg.p,Xi,X_ );
[ E_ ] = getSubDivKVValues( Eta, 1);
[B,Eta]=RefineKnotSolidEta( B,deg.q,Eta,E_ );
[ Z_ ] =getSubDivKVValues( Zeta, 1);
[B,Zeta]=RefineKnotSolidZeta( B,deg.r,Zeta,Z_ );

% Number of control points after refinement
n = size(B,1);
m = size(B,2);
l = size(B,3);

%% Build connectivity arrays
% Knot vectors in analysis
KV.Xi =Xi; KV.Eta = Eta; KV.Zeta = Zeta; clear Xi Eta Zeta

% Build connectivity arrays
nel = (n-deg.p) * (m-deg.q) * (l-deg.r); % number of elements
nnp = n*m*l; % number of global basis functions
nen = (deg.p+1)*(deg.q+1)*(deg.r+1); % number of local basis functions
ndof = nnp*3; % number of global degrees of freedom
ldof = nen*3; % number of local degrees of freedom
% Build connectivity arrays
[INN,IEN] = BldINCIEN( deg.p,deg.q,deg.r,n,m,l ); % = INC (is unique for IGA), IEN(=ENOD', if there is a CALFEM counterpart it might be called ENOD, but transposed)
ID = reshape(1:max(max(IEN))*3,3,max(max(IEN))); % ~= DOF' (Similar to DOF in CALFEM, but transposed)
LM = zeros(3*nen,nel); % ~= EDOF' (Similar to EDOF in CALFEM, but transposed)
for i = 1 : nel
    LM(:,i)=reshape(ID(:,IEN(:,i)),3*nen,1);
end

%% Material parameters:
E_Y = 210e9;
nu_P = 0.3;
D=hooke(4,E_Y,nu_P);

%% Gauss-Legendre quadrature points:
[ gp_x,w_x ] = getGP( deg.p );
[ gp_y,w_y ] = getGP( deg.q );
[ gp_z,w_z ] = getGP( deg.r );
NQUADx = size(gp_x,2);
NQUADy = size(gp_y,2);
NQUADz = size(gp_z,2);

%% Stiffness matrix and load vector computation

% Element loop
K = zeros(ndof); % Needs to be changed to sparse for large problems!!
% K = zeros(ndof);
F = zeros(ndof,1);
for e = 1 : nel
    % NURBS coordinates; convention consistent with Algorithm 7 in Hughes
    ni = INN(IEN(1,e),1);
    nj = INN(IEN(1,e),2);
    nk = INN(IEN(1,e),3);
    
    % Check if element has zero measure
    if (KV.Xi(ni+1) == KV.Xi(ni)) || (KV.Eta(nj+1) == KV.Eta(nj)) || (KV.Zeta(nk+1) == KV.Zeta(nk))
        continue
    end
    
    Ke = zeros(nen*3);
    Fe = zeros(nen*3,1);
    for i = 1 : NQUADx % Loop trough Gauss points
        for j = 1 : NQUADy
            for k = 1 : NQUADz
                % Gauss point
                GP.xi_tilde = gp_x(i);
                GP.eta_tilde = gp_y(j);
                GP.zeta_tilde = gp_z(k);
                
                % Get Basis, derivatives, and det(J) for current gauss pt
                [ ~,dR_dx,Jdet ] = Shape_function( GP,e,deg,B,KV,INN,IEN);
                
                % Combine quadrature weights with det(J)
                Jmod = abs(Jdet)*w_x(i)*w_y(j)*w_z(k);
                
                % Build Ke
                [ Ke_ ] = Build_K_Local( dR_dx,Jmod,D,nen );
                Ke = Ke + Ke_;
            end
        end
    end

    % Global Assembly
    idx = LM(:,e)';
    K(idx,idx) = K(idx,idx) + Ke;
    F(idx) = F(idx)+Fe;
end

% Apply load at two corner nodes
loc_z=0.1;
loc_x=3.0;
constNod = [];
for i = 1 : numel(B)
    if B{i}(3) == loc_z
        if (B{i}(2) == 0.0) || (B{i}(2) == 0.1)
            if B{i}(1) == loc_x
                constNod=[constNod i];
            end
        end
        
    end
end
% Apply load in vertical direction on identified nodes that are listed in
% constNod:
F(ID(3,constNod)) = -3e4;

%% Boundary condiitons
% Find controlpoints with x = 0 to constrain
loc_x=0;
constNod = [];
for i = 1 : numel(B)
    if B{i}(1) == loc_x
        constNod=[constNod i];
    end
end
bc1=reshape(ID(:,constNod),numel(ID(:,constNod)),1);
bc=[bc1,zeros(length(bc1),1)];

%% Solve system
a=solveq(K,F,bc);

% Rewrite to cell structure. (For plotting)
u = cell(size(B));
comb = cell(size(B));
for i = 1 : size(ID,2)
    u{i} = [a(ID(:,i)); 0];
    comb{i} = B{i} + u{i};
end

%% plot
plotNurbsSolidElementSimple( KV, B )
title('Initial geometry')
plotNurbsSolidElementSimple( KV, comb )
title('Displacements')




