function [  ] = plotNurbsSolidElement( Cx, Cy, Cz, Xi, Eta, Zeta, U, V, W )
%[  ] = plotNurbsSolidElement( Cx, Cy, Cz, Xi, Eta, Zeta, U, V, W )
%-------------------------------------------------------------
% PURPOSE:
%   Plot NURBS solid geometry, and draw knots (elements).
%
%
% INPUT: Cx = X-coordinates at each evaluation point of basis functions (size : (n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1),(n-r+1 x resz*nnzKz+1))
%        Cy = Y-coordinates at each evaluation point of basis functions (size : (n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1),(n-r+1 x resz*nnzKz+1))
%        Cz = Z-coordinates at each evaluation point of basis functions (size : (n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1),(n-r+1 x resz*nnzKz+1))
%        Xi = non-uniform clamped knot vector
%        Eta = non-uniform clamped knot vector
%        Zeta = non-uniform clamped knot vector
%        U = corresponding parameter values to Cx
%        V = corresponding parameter values to Cy
%        W = corresponding parameter values to Cz
%
% OUTPUT: none
%-------------------------------------------------------------

title('Elements')
view([60 10])
grid on
% Front:
surfl(Cx(:,:,1),Cy(:,:,1),Cz(:,:,1))
hold on
% Back:
surfl(Cx(:,:,end),Cy(:,:,end),Cz(:,:,end))

% Side1:
Ax(:,:) = Cx(:,end,:);
Ay(:,:) = Cy(:,end,:);
Az(:,:) = Cz(:,end,:);
surfl(Ax,Ay,Az)

% Side 2:
Ax(:,:) = Cx(:,1,:);
Ay(:,:) = Cy(:,1,:);
Az(:,:) = Cz(:,1,:);
surfl(Ax,Ay,Az)

% Bottom:
clear Ax Ay Az
Ax(:,:) = Cx(1,:,:);
Ay(:,:) = Cy(1,:,:);
Az(:,:) = Cz(1,:,:);
surfl(Ax,Ay,Az)

% Top:
clear Ax Ay Az
Ax(:,:) = Cx(end,:,:);
Ay(:,:) = Cy(end,:,:);
Az(:,:) = Cz(end,:,:);
surfl(Ax,Ay,Az)

shading interp
colormap summer
hidden off
axis equal

%% Nurbs surface basis from given knots and weigths(in Control point data).

% Find knot positions
Xiu=unique(Xi);
xi_pos=[];
for i = 1 : length(Xiu)
    xi_pos = [xi_pos find(U == Xiu(i))];
end
% Find knot positions
Etau=unique(Eta);
eta_pos=[];
for i = 1 : length(Etau)
    eta_pos = [eta_pos find(V == Etau(i))];
end
% Find knot positions
Zetau=unique(Zeta);
zeta_pos=[];
for i = 1 : length(Zetau)
    zeta_pos = [zeta_pos find(W == Zetau(i))];
end

col = [0.5020    0.7490    0.5529];
% Loop trough knot positions and plot
for i = xi_pos
    for j = eta_pos
        for k = zeta_pos
            clear x y z
            x(:,1) = Cx(i,j,:);
            y(:,1) = Cy(i,j,:);
            z(:,1) = Cz(i,j,:);
            plot3(x,y,z,'LineWidth',2,'Color',col)
            clear x y z
            x(:,1) = Cx(i,:,k);
            y(:,1) = Cy(i,:,k);
            z(:,1) = Cz(i,:,k);
            plot3(x,y,z,'LineWidth',2,'Color',col)
            clear x y z
            x(:,1) = Cx(:,j,k);
            y(:,1) = Cy(:,j,k);
            z(:,1) = Cz(:,j,k);
            plot3(x,y,z,'LineWidth',2,'Color',col)
        end
    end
end

end

