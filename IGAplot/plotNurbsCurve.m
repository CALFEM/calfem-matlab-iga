function [  ] = plotNurbsCurve( Cx, Cy, Cz, Xi, U, x, y, z )
%[  ] = plotNurbsCurve2D( Cx, Cy, Xi, U, x, y )
%-------------------------------------------------------------
% PURPOSE:
% Plots the NURBS basis curve in 3D in R, U. To get R, U use
% the function nrbasis_num
%
% INPUT: Cx = x-coordinates (1 x res*nnzK+1)
%        Cy = y-coordinates (1 x res*nnzK+1)
%        Cz = z-coordinates (1 x res*nnzK+1)
%        U = corresponding KV parameter values (1 x res*nnzK+1)
%        Xi = Knot vector
%        x = x-coordinates of control points
%        y = y-coordinates of control points
%        z = z-coordinates of control points
%
% OUTPUT: none
%-------------------------------------------------------------


plot3(Cx,Cy,Cz,'color',[0 0.4470 0.7410]);
grid on
hold on
plot3(x,y,z,'r:o')

%% plot knots
Xiu = unique(Xi);
xi_pos=[];
for i = 1 : length(Xiu)
    xi_pos = [xi_pos find(U == Xiu(i))];
end
plot3(Cx(xi_pos),Cy(xi_pos),Cz(xi_pos),'d','color',[0 0.4470 0.7410])
axis equal


legend('Nurbs curve','Control points','Knots','Location','Best')

end

