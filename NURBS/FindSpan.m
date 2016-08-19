function knotSpanIndex = FindSpan(n,p,u,U)
% knotSpanIndex = FindSpan(n,p,u,U)
%-------------------------------------------------------------
% PURPOSE:
%  find the knot span index for one variable u, NURBS-Book 
%  (algorithm A2.1)
% 
% INPUT: n = number of basis function - 1
%        p = degree of the basis functions
%        u = evaluation point
%        U = knot vector (row vector)
%
% OUTPUT: knotSpanIndex = index of knot span
%-------------------------------------------------------------

if (u == U(n+2))
    knotSpanIndex= n-1;
    return
end
low = p;
high = n+1;
mid = floor((low + high)/2);
while (u <U(mid+1) || u >= U(mid+2) )
    if( u < U(mid+1))
        high = mid;
    else
        low = mid;
    end
    mid = floor((low+high)/2);
end
knotSpanIndex = mid;
end