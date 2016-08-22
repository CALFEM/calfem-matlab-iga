function [  ] = plotNurbsSolidElementSimple( KV, B )
%[  ] = plotNurbsSolidElementSimple( KV, B )
%-------------------------------------------------------------
% PURPOSE:
%   Plot NURBS solid geometry, and draw knots (elements).
%
%
% INPUT: KV = Knot vector with:
%        KV.Xi, KV.Eta, KV.Zeta
%        B = Control point net (Cell)
%
% OUTPUT: none
%-------------------------------------------------------------

Xi = KV.Xi; Eta = KV.Eta; Zeta = KV.Zeta;
clear KV.Xi KV.Eta KV.Zeta

if numel(B) < 300
res_x = 5;
res_y = 5;
res_z = 5;
else
res_x = 1;
res_y = 1;
res_z = 1;    
end

% Nurbs surface basis from given knots and weigths(in Control point data).
[R, U, V, W] = nrbasis_solid_num(Xi,Eta,Zeta,B,res_x,res_y,res_z);

% Deformed NURBS solid according to (2.29) in Cotrell, Hughes & Bazilevs
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

%% Plotting the deformed Nurbs solid
figure()
plotNurbsSolidElement( Cx, Cy, Cz, Xi, Eta, Zeta, U, V, W )

end

