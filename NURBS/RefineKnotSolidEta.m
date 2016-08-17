function [ Q,Etai ] = RefineKnotSolidEta( B,q,Eta,E )
%[ Q,Etai ] = RefineKnotSolidEta( B,q,Eta,E )
%-------------------------------------------------------------
% PURPOSE:
% Performs knot refinement for NURBS surfaces and solids in the 
% second direction (Eta).
%
% INPUT: B = cell structure containing 3-D control points and weigths
%            B{i,j,k}[x,y,z,w]
%
%        q = degree
%        Eta = non-uniform clamped knot vector
%        E = Knots to instert into Xi. 
%
% OUTPUT: Q = new cell structure containing 3-D control points 
%             and weights
%         Etai = new knot vector
%-------------------------------------------------------------

Q = cell([size(B,1),size(B,2)+length(E),size(B,3)]);

for i = 1 : size(B,1)
    for k = 1 : size(B,3)
        for j = 1 : size(B,2)
            Pw(j,:) = [B{i,j,k}(1)*B{i,j,k}(4) B{i,j,k}(2)*B{i,j,k}(4) B{i,j,k}(3)*B{i,j,k}(4) B{i,j,k}(4) ];
        end
        r_ = length(E)-1;
        [Eta_bar,Qw] = RefineKnotVectCurve(size(Pw,1)-1,q,Eta,Pw,E,r_);
        qw = Qw(:,4);
        Qw(:,1) = Qw(:,1)./qw;
        Qw(:,2) = Qw(:,2)./qw;
        Qw(:,3) = Qw(:,3)./qw;
        for j = 1 : size(B,2)+length(E)
            Q{i,j,k}(1:4,1) = Qw(j,1:4);
        end
    end
end

Etai = Eta_bar;
end

