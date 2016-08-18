function [  ] = plotNurbsSolid( Cx,Cy,Cz, B )
%[  ] = plotNurbsSolid( Cx,Cy,Cz, B )
%-------------------------------------------------------------
% PURPOSE:
%   Plot NURBS solid geometry.
%
%
% INPUT: Cx = X-coordinates at each evaluation point of basis functions (size : (n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1),(n-r+1 x resz*nnzKz+1))
%        Cy = Y-coordinates at each evaluation point of basis functions (size : (n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1),(n-r+1 x resz*nnzKz+1))
%        Cz = Z-coordinates at each evaluation point of basis functions (size : (n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1),(n-r+1 x resz*nnzKz+1))
%        B = Cell structure containing control points
%
% OUTPUT: none
%-------------------------------------------------------------



%Front:
surf(Cx(:,:,1),Cy(:,:,1),Cz(:,:,1))
grid on
hold on

%Back:
surf(Cx(:,:,end),Cy(:,:,end),Cz(:,:,end))

% Side 1:
Ax(:,:) = Cx(:,end,:);
Ay(:,:) = Cy(:,end,:);
Az(:,:) = Cz(:,end,:);
surf(Ax,Ay,Az)

% Side 2:
Ax(:,:) = Cx(:,1,:);
Ay(:,:) = Cy(:,1,:);
Az(:,:) = Cz(:,1,:);
surf(Ax,Ay,Az)

% Bottom:
clear Ax Ay Az
Ax(:,:) = Cx(1,:,:);
Ay(:,:) = Cy(1,:,:);
Az(:,:) = Cz(1,:,:);
surf(Ax,Ay,Az)

% Top:
clear Ax Ay Az
Ax(:,:) = Cx(end,:,:);
Ay(:,:) = Cy(end,:,:);
Az(:,:) = Cz(end,:,:);
surf(Ax,Ay,Az)

xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%% Plotting the control polygon
% Build points
% for i = 1: size(B,1)
%     for j = 1 : size(B,2);
%         for k = 1 : size(B,3);
%             x(i,j,k) = B{i,j,k}(1);
%             y(i,j,k) = B{i,j,k}(2);
%             z(i,j,k) = B{i,j,k}(3);
%         end
%     end
% end
% for i = 1 : size(x,3)
% plot3(x(:,:,i),y(:,:,i),z(:,:,i),'o')
% mesh(x(:,:,i),y(:,:,i),z(:,:,i))
% end
% for i = 1 : size(x,2)
%     x_(:,:)=x(:,i,:);
%     y_(:,:)=y(:,i,:);
%     z_(:,:)=z(:,i,:);
%     mesh(x_,y_,z_)
% end
% clear x_ y_ z_
% for i = 1 : size(x,1)
%     x_(:,:)=x(i,:,:);
%     y_(:,:)=y(i,:,:);
%     z_(:,:)=z(i,:,:);
%     mesh(x_,y_,z_)
% end
% hidden off
% axis equal
% title('Nurbs solid')

end

