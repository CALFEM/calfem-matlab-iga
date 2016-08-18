function [ X_ ] = getIntermid( Xi )
% [ X_ ] = getIntermid( Xi )
%-------------------------------------------------------------
% PURPOSE
%  Computes intermediate values in a knot vector for a single sub-division 
%
% INPUT: Xi : Univariate Knot-Vector
%
% OUTPUT:  X_ : Parameter values for sub-division
%-------------------------------------------------------------

X_ = zeros(1,length(Xi)-1);
for i = 1 : numel(X_)
    X_(i) = (Xi(i) + Xi(i+1))/2;
end
X_ = setdiff(X_,Xi);
end

