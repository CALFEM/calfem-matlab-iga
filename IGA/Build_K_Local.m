function [ Ke ] = Build_K_Local( dR_dx,Jmod,D,nen )
% [ Ke ] = Build_K_Local( dR_dx,Jmod,D,nen )
%-------------------------------------------------------------
% PURPOSE:
% Calculate the contribution from the current integration point 
% to the linear elastic stiffness matrix for a solid 
% NURBS element.
%
% INPUT: dR_dx = vector of basis function derivatives (nen x 3)
%
%        Jmod  = Jacobian determinant
%
%        D     = constitutive matrix (6 x 6)
%
%        nen   = number of local basis functions
%
% OUTPUT: Ke = element stiffness matrix contribution (nen*3 x nen*3)
%-------------------------------------------------------------

% Generate B matrix
B=zeros(6,nen*3);
B(1,1:3:end) = dR_dx(:,1);
B(2,2:3:end) = dR_dx(:,2);
B(3,3:3:end) = dR_dx(:,3);
B(4,2:3:end) = dR_dx(:,3);
B(4,3:3:end) = dR_dx(:,2);
B(5,1:3:end) = dR_dx(:,3);
B(5,3:3:end) = dR_dx(:,1);
B(6,1:3:end) = dR_dx(:,2);
B(6,2:3:end) = dR_dx(:,1);

% Contrubution to element stiffness matrix
Ke = B'*D*B*Jmod;

end

