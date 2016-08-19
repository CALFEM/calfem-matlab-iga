function [ Q,Zetai ] = DegElevZeta( B,r,Zeta,d )
%[ Q,Zetai ] = DegElevZeta( B,r,Zeta,d )
%-------------------------------------------------------------
% PURPOSE:
%  Performs degree elevation for NURBS solids in the 
%  third direction (Zeta).
%
% INPUT: B = cell structure containing 3-D control points and weigths
%            B{i,j,k}[x,y,z,w]
%
%        r = degree
%        Zeta = non-uniform clamped knot vector
%        d = number of times to raise degree
%
% OUTPUT: Q = new cell structure containing 3-D control points 
%             and weights
%         Zetai = new knot vector
%-------------------------------------------------------------

for i = 1 : size(B,1)
    for j = 1 : size(B,2)
        for k = 1 : size(B,3)
            Pw(k,:) = [B{i,j,k}(1)*B{i,j,k}(4) B{i,j,k}(2)*B{i,j,k}(4) B{i,j,k}(3)*B{i,j,k}(4) B{i,j,k}(4) ];
        end
        [Qw,Zeta_bar] = bspdegelev(r,Pw',Zeta,d);
        Qw=Qw';
        qw = Qw(:,4);
        Qw(:,1) = Qw(:,1)./qw;
        Qw(:,2) = Qw(:,2)./qw;
        Qw(:,3) = Qw(:,3)./qw;
        for k = 1 : size(Qw,1)
            Q{i,j,k}(1:4,1) = Qw(k,1:4);
        end
    end
end

Zetai = Zeta_bar;
end

