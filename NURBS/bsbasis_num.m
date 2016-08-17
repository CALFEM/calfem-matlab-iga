function [ N, Xi_store ] = bsbasis_num( Xi, p, res )
%[ NN_, U_ ] = bsbasis_num( Xi, p, res )
%-------------------------------------------------------------
% PURPOSE:
% Returns all B-Spline basis functions of order p for the k 
% intervals in the knot vector Xi, and n weghts in w.
%
% INPUT: Xi = non-uniform clamped knot vector (1 x k+1)
%        p = degree
%        res = number of parameter values per non-zero knot span (nnzK)
%
% OUTPUT: N = basis functions (n-p+1 x res*nnzK+1)
%         Xi_store = corresponding parameter values (1 x res*nnzK+1)
%-------------------------------------------------------------

k = length(Xi);
n = k - p - 1;

Xi_store =  min(Xi) : max(Xi)/res : max(Xi); % Evaluate basis functions at each xi in Xi_store
Xi_store = unique(sort([Xi_store Xi])); % <- Makes sure that the knots are included (In order to plot the knots, otherwise we could ignore this)
N = zeros(n,length(Xi_store)); % We have n basisfunctions, evaluated at each xi
for i = 1:length(Xi_store)
    xi = Xi_store(i); % Get current xi
    j=FindSpan(n,p,xi,Xi); % Find knot span
    N_loc=BasisFun(j,xi,p,Xi); % Get only the non-zero basis functions (there are p + 1 non-zero basis functions).
    N(j-p+1:j+1,i)=N_loc; % Store the non-zero basis functions
end

end