function [  ] = plotNurbsBasis( R, U )
%[  ] = plotNurbsBasis( R, U )
%-------------------------------------------------------------
% PURPOSE:
% Plots the NURBS basis functions in R, U. To get R, U use
% the function nrbasis_num
%
% INPUT: R = basis functions (n-p+1 x res*nnzK+1)
%        U = corresponding KV parameter values (1 x res*nnzK+1)
%
% OUTPUT: none
%-------------------------------------------------------------

c = char('kmcrgb');
l = 1;
for j = 1 : size(R,1)
    plot(U,R(j,:),c(l));
    hold on
    
    l = l +1;
    if l > length(c)
        l = 1;
    end
end
grid on
end

