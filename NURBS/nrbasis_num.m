function [ R , U] = nrbasis_num( Xi, w,res )
%[ R , U] = nrbasis_num( Xi, w,res )
%-------------------------------------------------------------
% PURPOSE:
%  Returns all NURBS basis functions of order p for the k 
%  intervals in the knot vector Xi, and n weghts in w.
%
% INPUT: Xi = non-uniform clamped knot vector (1 x k+1)
%        w = weights of controlpoints (1 x n)
%        res = number of parameter samples in total knot vector, 
%              uniformly distributed (nnzK)
%
% OUTPUT: R = basis functions (n-p+1 x res*nnzK+1)
%         U = corresponding parameter values (1 x res*nnzK+1)
%-------------------------------------------------------------

% Number of knots
k = length(Xi);

% Number of weights
n = length(w);

% Order of basis
p = k-n-1;

% B-spline basis, returns N and correspding parameter values U
[ N, U ] = bsbasis_num( Xi, p, res );

% Weight W(xi) according to (2.27) in Cotrell, Hughes & Bazilevs
W = N'*w;

% Nurbs Basis according to (2.27) in Cotrell, Hughes & Bazilevs
R = zeros(size(N));
for i = 1 : size(N,2)
    R(:,i) = N(:,i).*w(:)./W(i);
end


end

