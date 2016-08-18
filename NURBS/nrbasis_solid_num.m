function [ R, U, V, W_ ] = nrbasis_solid_num( Xi, Eta, Zeta, B, resx, resy, resz )
%[ R, U, V, W_ ] = nrbasis_solid_num( Xi, Eta, Zeta, B, resx, resy, resz ))
%-------------------------------------------------------------
% PURPOSE:
%   Returns the basis functions R(xi) of order p,q for the a intervals in
%   xi_i = [xi_1, xi_2, ... x_a+1],  and b intervals in
%   eta_i =[eta_1, ..., eta_b+1],    and c intervals in
%   zeta_i =[zeta_1, ..., zeta_c+1], and m x n x l weights in
%   B{m,n,l} = [x,y,z,w_mnl]
%
%
% INPUT: Xi = non-uniform clamped knot vector (1 x a+1)
%        Eta = non-uniform clamped knot vector (1 x b+1)
%        Zeta = non-uniform clamped knot vector (1 x c+1)
%        B{m,n,l} = [x,y,z,w_mn]
%        resx = number of parameter values per non-zero knot span (nnzKx)
%        resy = number of parameter values per non-zero knot span (nnzKy)
%        resz = number of parameter values per non-zero knot span (nnzKz)
%
% OUTPUT: R{m,n,l} = basis functions {(n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1),(n-r+1 x resz*nnzKz+1)}
%         U = corresponding parameter values (1 x resx*nnzKx+1)
%         V = corresponding parameter values (1 x resy*nnzKy+1)
%         W = corresponding parameter values (1 x resz*nnzKz+1)
%-------------------------------------------------------------

% Number of knots
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);

% Number of weights
n = size(B,1);
m = size(B,2);
l = size(B,3);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;

% B-spline basis
[ N_ip, U ] = bsbasis_num( Xi, p, resx );
[ N_jq, V ] = bsbasis_num( Eta, q, resy );
[ N_kr, W_ ] = bsbasis_num( Zeta, r, resz );

% Weight W(xi) according to (2.27) in Cotrell, Hughes & Bazilevs
W = zeros(size(N_ip,2),size(N_jq,2),size(N_kr,2));
for i = 1 : n
    for j = 1 : m  
        for k = 1 : l
            X = N_ip(i,:).'*N_jq(j,:);
            Y = reshape(X(:)*N_kr(k,:), [ size(X,1) size(X,2) length(N_kr(k,:)) ]);
            W = W + Y*B{i,j,k}(4);
        end
    end
end

% Nurbs Basis according to (2.27) in Cotrell, Hughes & Bazilevs
R = cell(i,j,k);
for i = 1 : size(N_ip,1)
    for j = 1 : size(N_jq,1)
        for k = 1 : size(N_kr,1)
            X = N_ip(i,:)'*N_jq(j,:);
            Y = reshape(X(:)*N_kr(k,:), [ size(X,1) size(X,2) length(N_kr(k,:)) ]);
            R{i,j,k} = Y*B{i,j,k}(4)./W;
        end
    end
end
end