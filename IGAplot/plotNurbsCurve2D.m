function [  ] = plotNurbsCurve2D( Cx, Cy, Xi, U, x, y )
%[  ] = plotNurbsCurve2D( Cx, Cy, Xi, U, x, y )
%-------------------------------------------------------------
% PURPOSE:
% Plots the NURBS basis curve in 2D in R, U. To get R, U use
% the function nrbasis_num
%
% INPUT: Cx = x-coordinates (1 x res*nnzK+1)
%        Cy = y-coordinates (1 x res*nnzK+1)
%        U = corresponding KV parameter values (1 x res*nnzK+1)
%        Xi = Knot vector
%        x = x-coordinates of control points
%        y = y-coordinates of control points
%
% OUTPUT: none
%-------------------------------------------------------------


plot(Cx,Cy,'color',[0 0.4470 0.7410]);
grid on
hold on
plot(x,y,'r:o')

%% plot knots
Xiu = unique(Xi);
xi_pos=[];
for i = 1 : length(Xiu)
    xi_pos = [xi_pos find(U == Xiu(i))];
end
plot(Cx(xi_pos),Cy(xi_pos),'d','color',[0 0.4470 0.7410])
axis equal


legend('Nurbs curve','Control points','Knots','Location','Best')

end

