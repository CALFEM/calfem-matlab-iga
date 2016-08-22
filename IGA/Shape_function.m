function [ R,dR_dx,J ] = Shape_function( GP,e,deg,B,KV,INC,IEN)
% [ R,dR_dx,J ] = Shape_function( GP,e,deg,B,KV,INN,IEN)
%-------------------------------------------------------------
% PURPOSE:
% For a solid NURBS element, calculate the the vector of local 
% shape functions R and an array of their derivatives dR_dx at
% the current gauss point, and the Jacobian determinant J.
% 
% INPUT: The quadrature point GP, containing:
%        GP.Xi_tilde, GP.Eta_tilde, and GP.Zeta_tilde 
%
%        element number e, 
%
%        polynomial orders deg, containing:
%        deg.p, deg.q, and deg.r
%
%        control net B{n,m,l}(x,y,z,w), with a cell structure 
%        with n points in Xi direction, m points in Eta 
%        direction, and l points in Zeta direction. Weights 
%        are stored as the fourth component.
%
%        knot vectors KV, containing the clamped non-uniform 
%        knot vectors: 
%        KV.Xi, KV.Eta, and KV.Zeta
%
%        and connectivity arrays INC and IEN.
% 
% OUTPUT: R = vector of basis functions (nen x 1)
%         dR_dx = vector of basis function derivatives (nen x 3)
%         J = Jacaboian determinant (1)
%-------------------------------------------------------------

p = deg.p; q = deg.q; r = deg.r;

% number of local basis functions:
nen = (p+1)*(q+1)*(r+1); 
                          
% NURBS coordinates; convention consistent with Algorithm 7
ni = INC(IEN(1,e),1);
nj = INC(IEN(1,e),2);
nk = INC(IEN(1,e),3);

% Calculate parametric coordinates from parent element coordinates
% Knot vectors KV.Xi, KV.Eta, and KV.Zeta and
% parent element coordinates xi_tilde, eta_tilde, zeta_tilde
% are given as input
xi = ((KV.Xi(ni+1)-KV.Xi(ni))*GP.xi_tilde ...
+ (KV.Xi(ni+1)+KV.Xi(ni))) / 2;
eta = ((KV.Eta(nj+1)-KV.Eta(nj))*GP.eta_tilde ...
+ (KV.Eta(nj+1)+KV.Eta(nj))) / 2;
zeta = ((KV.Zeta(nk+1)-KV.Zeta(nk))*GP.zeta_tilde ...
+ (KV.Zeta(nk+1)+KV.Zeta(nk))) / 2;


% Calculate univariate B-spline functions using (2.1) and (2.2)
% and their derivatives using (2.12)
N1 = Der1BasisFun(ni-1,xi,p,KV.Xi)'; % xi-dir.
N2 = Der1BasisFun(nj-1,eta,q,KV.Eta)'; % eta-dir.
N3 = Der1BasisFun(nk-1,zeta,r,KV.Zeta)'; % zeta-dir.
N = N1(:,1);
dN_dxi = N1(:,2);
M = N2(:,1);
dM_deta = N2(:,2);
L = N3(:,1);
dL_dzeta = N3(:,2);
clear N1 N2 N3

% Build numerators and denominators (in local numbering)
x = zeros(1,nen); y = zeros(1,nen); z = zeros(1,nen);
R = zeros(nen,1);   % Array of trivariate NURBS basis functions
dR_dxi = zeros(nen,3); % Trivariate NURBS function derivatives
                       % w.r.t. parametric coordinates
loc_num = 0; % Local basis function counter
for k = 0 : r
    for j = 0 : q
        for i = 0 : p
            loc_num = loc_num + 1;
            
            R(loc_num) = N(p+1-i)*M(q+1-j)*L(r+1-k) ...
                         * B{ni-i,nj-j,nk-k}(4); % Function numerator (N*M*L*w)
            
            % Get coordinates in local numbering
            x(loc_num) = B{ni-i,nj-j,nk-k}(1);
            y(loc_num) = B{ni-i,nj-j,nk-k}(2);
            z(loc_num) = B{ni-i,nj-j,nk-k}(3);
            %w(loc_num) = B{ni-i,nj-j,nk-k}(4);
            
            dR_dxi(loc_num,1) = dN_dxi(p+1-i)*M(q+1-j)*L(r+1-k) ...
            * B{ni-i,nj-j,nk-k}(4); % Derivative numerator (dN*M*L*w)
            
            dR_dxi(loc_num,2) = N(p+1-i)*dM_deta(q+1-j)*L(r+1-k) ...
            * B{ni-i,nj-j,nk-k}(4); % Derivative numerator (N*dM*L*w)

            dR_dxi(loc_num,3) = N(p+1-i)*M(q+1-j)*dL_dzeta(r+1-k) ...
            * B{ni-i,nj-j,nk-k}(4); % Derivative numerator (N*M*dL*w)
        end
    end
end
W = sum(R); % Function denominator (Sum(N*M*L*w))
dW_dxi = sum(dR_dxi(:,1)); % Derivative denominator (Sum(dN*M*L*w))
dW_deta = sum(dR_dxi(:,2)); % Derivative denominator (Sum(N*dM*L*w))
dW_dzeta = sum(dR_dxi(:,3)); % Derivative denominator (Sum(N*M*dL*w))
            
% Divide by denominators to complete definitions of functions
% and derivatives w.r.t. parametric coordinates
dR_dxi(:,1) = (dR_dxi(:,1)*W - R*dW_dxi) / W^2;
dR_dxi(:,2) = (dR_dxi(:,2)*W - R*dW_deta) / W^2;
dR_dxi(:,3) = (dR_dxi(:,3)*W - R*dW_dzeta) / W^2;
R = R/W;

% Gradient of mapping from parameter space to physical space
dx_dxi = [x;y;z] * dR_dxi;

% Compute derivatives of basis functions
% with respect to physical coordinates
dR_dx = dR_dxi/dx_dxi; %dR_dxi * inv(dx_dxi)

% Gradient of mapping from parent element to parameter space
dxi_dtildexi=zeros(3); % Derivative of parametric coordinates
                       % w.r.t. parent element coordinates
dxi_dtildexi(1,1) = (KV.Xi(ni+1)-KV.Xi(ni))/2;
dxi_dtildexi(2,2) = (KV.Eta(nj+1)-KV.Eta(nj))/2;
dxi_dtildexi(3,3) = (KV.Zeta(nk+1)-KV.Zeta(nk))/2;

% Compute the Jacobian
J_mat = dx_dxi*dxi_dtildexi;

% Compute Jacobian determinant
J = det(J_mat);

end

