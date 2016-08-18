function [ R, U, V ] = nrbasis_surf_num( Xi, Eta, B, resx, resy )
%[ R, U, V ] = nrbasis_surf_num( Xi, Eta, B, resx, resy )
%-------------------------------------------------------------
% PURPOSE:
%   Returns the basis functions R(xi) of order p,q for the k intervals in
%   xi_i = [xi_1, xi_2, ... x_k+1], and l intervals in
%   eta_i =[eta_1, ..., eta_l+1],   and m x n weights in
%   B{m,n} = [x,y,z,w_mn]
%
%
% INPUT: Xi = non-uniform clamped knot vector (1 x k+1)
%        Eta = non-uniform clamped knot vector (1 x l+1)
%        B{m,n} = [x,y,z,w_mn]
%        resx = number of parameter values per non-zero knot span (nnzKx)
%        resy = number of parameter values per non-zero knot span (nnzKy)
%
% OUTPUT: R{m,n} = basis functions {(n-p+1 x resx*nnzKx+1),(m-q+1 x resy*nnzKy+1)}
%         U = corresponding parameter values (1 x resx*nnzKx+1)
%         V = corresponding parameter values (1 x resy*nnzKy+1)
%-------------------------------------------------------------

% Number of knots
k = length(Xi);
l = length(Eta);

% Number of weights
m = size(B,1);
n = size(B,2);

% Order of basis
p = k-m-1;
q = l-n-1;

% B-spline basis
[ N_ip, U ] = bsbasis_num( Xi, p, resx );
[ N_jq, V ] = bsbasis_num( Eta, q, resy );

% Weight W(xi) according to (2.27) in Cotrell, Hughes & Bazilevs
W = zeros(size(N_ip,2),size(N_jq,2));
for i = 1 : m
    for j = 1 : n  
        W = W + N_ip(i,:).'*N_jq(j,:)*B{i,j}(4);
    end
end

% Nurbs Basis according to (2.27) in Cotrell, Hughes & Bazilevs
R = cell(i,j);
for i = 1 : size(N_ip,1)
    for j = 1 : size(N_jq,1)
        R{i,j} = N_ip(i,:)'*N_jq(j,:)*B{i,j}(4)./W;
        R{i,j} = sparse(R{i,j});
    end
end
end