function [ Q,Xii ] = RefineKnotSolidXi( B,p,Xi,X )
%[ Q,Xii ] = RefineKnotSolidXi( B,p,Xi,X )
%-------------------------------------------------------------
% PURPOSE:
% Performs knot refinement for NURBS surfaces and solids in the 
% first direction (Xi).
%
% INPUT: B = cell structure containing 3-D control points and weigths
%            B{i,j,k}[x,y,z,w]
%
%        p = degree
%        Xi = non-uniform clamped knot vector
%        X = Knots to instert into Xi. 
%
% OUTPUT: Q = new cell structure containing 3-D control points 
%             and weights
%         Xii = new knot vector
%-------------------------------------------------------------

Q = cell([size(B,1)+length(X),size(B,2),size(B,3)]);

for j = 1 : size(B,2)
    for k = 1 : size(B,3)
        for i = 1 : size(B,1)
            Pw(i,:) = [B{i,j,k}(1)*B{i,j,k}(4) B{i,j,k}(2)*B{i,j,k}(4) B{i,j,k}(3)*B{i,j,k}(4) B{i,j,k}(4) ];
        end
        r_ = length(X)-1;
        [Xi_bar,Qw] = RefineKnotVectCurve(size(Pw,1)-1,p,Xi,Pw,X,r_);
        qw = Qw(:,4);
        Qw(:,1) = Qw(:,1)./qw;
        Qw(:,2) = Qw(:,2)./qw;
        Qw(:,3) = Qw(:,3)./qw;
        for i = 1 : size(B,1)+length(X)
            Q{i,j,k}(1:4,1) = Qw(i,1:4);
        end
    end
end

Xii = Xi_bar;
end

