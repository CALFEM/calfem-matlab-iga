function [ Q,Zetai ] = RefineKnotSolidZeta( B,r,Zeta,Z )
%[ Q,Zetai ] = RefineKnotSolidZeta( B,r,Zeta,Z )
%-------------------------------------------------------------
% PURPOSE:
% Performs knot refinement for NURBS solids in the 
% third direction (Zeta).
%
% INPUT: B = cell structure containing 3-D control points and weigths
%            B{i,j,k}[x,y,z,w]
%
%        r = degree
%        Zeta = non-uniform clamped knot vector
%        Z = Knots to instert into Xi. 
%
% OUTPUT: Q = new cell structure containing 3-D control points 
%             and weights
%         Zetai = new knot vector
%-------------------------------------------------------------

Q = cell([size(B,1),size(B,2),size(B,3)+length(Z)]);

for i = 1 : size(B,1)
    for j = 1 : size(B,2)
        for k = 1 : size(B,3)
            Pw(k,:) = [B{i,j,k}(1)*B{i,j,k}(4) B{i,j,k}(2)*B{i,j,k}(4) B{i,j,k}(3)*B{i,j,k}(4) B{i,j,k}(4) ];
        end
        r_ = length(Z)-1;
        [Zeta_bar,Qw] = RefineKnotVectCurve(size(Pw,1)-1,r,Zeta,Pw,Z,r_);
        qw = Qw(:,4);
        Qw(:,1) = Qw(:,1)./qw;
        Qw(:,2) = Qw(:,2)./qw;
        Qw(:,3) = Qw(:,3)./qw;
        for k = 1 : size(B,3)+length(Z)
            Q{i,j,k}(1:4,1) = Qw(k,1:4);
        end
    end
end

Zetai = Zeta_bar;
end

