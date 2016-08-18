function [ X_div ] = getSubDivKVValues( Xi, sDivs )
% [ X_old ] = getSubDivKVValues( Xi, sDivs )
%-------------------------------------------------------------
% PURPOSE
%  Computes knot insertion points for sub-division
%
% INPUT: Xi : Univariate Knot-Vector
%        sDivs : number of sub-divisions
%
% OUTPUT:  X_div : Parameter values for sub-division
%-------------------------------------------------------------

[ X_ ] = getIntermid( Xi );
X_div = X_;
while sDivs > 1
    [ X_ ] = getIntermid( sort([Xi X_div]) );
    X_div = [X_div X_];
    sDivs = sDivs - 1;
end
X_div = sort(X_div);

end

