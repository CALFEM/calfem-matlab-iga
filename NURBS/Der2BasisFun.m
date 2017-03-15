function ders = Der2BasisFun(i,u,p,U)
%ders = Der2BasisFun(i,u,p,U)
%--------------------------------------------------------------
% PURPOSE:
% NURBS-Book modified (algorithm A2.3)
% evalute nonzero basis functions and first derivative
%
% INPUT:
%        i = current knotspan
%        u = evaluation point
%        p = degree of the basis functions
%        U = knot vector
% OUTPUT:
%        ders = matrix (2 , p+1)
%               1. row vector (dim p+1) values of the basis function N_(i-p) ... N_(i) at the evaluation point u
%               2. first derivation
%               3. second derivation
%--------------------------------------------------------------

ders   = zeros(2,p+1);
N      = zeros(p+1,p+1);
N(1,1) = 1;
left   = zeros(1,p+1);
right  = zeros(1,p+1);

for j=1:p
    left(j+1)  = u-U(i+1-j+1);
    right(j+1) = U(i+j+1)-u;
    saved = 0;
    for r=0:j-1
        N(j+1,r+1) = right(r+2) + left(j-r+1);
        temp = N(r+1,j)/N(j+1,r+1);
        N(r+1,j+1) = saved + right(r+2)*temp;
        saved = left(j-r+1)*temp;
    end
    N(j+1,j+1) = saved;
end

for j=0:p % Load the basis functions
    ders(1,j+1) = N(j+1,p+1);
end


%derivatives (k = 1 to 2)
for r=0:p %Loop over function index
    s1 = 0; s2 = 1;
    a(1,1) = 1;
    %kth derivative
    for k = 1 : 2
        d = 0.0;
        rk = r - k;
        pk = p -k;

        if (r >= k)
            a(s2+1,1) = a(s1+1,1)/N(pk+1+1,rk+1);
            d = a(s2+1,1) * N(rk+1,pk+1);
        end

        if (rk>=-1)
            j1=1;
        else
            j1 = -rk;
        end

        if(r-1 <= pk)
            j2 = k-1;
        else
            j2 = p-r;
        end
        
        for j = j1:j2
            a(s2+1,j+1) = ( a(s1+1,j+1)-a(s1+1,j) )/N(pk+1+1,rk+j+1);
            d = d + a(s2+1,j+1)*N(rk+j+1,pk+1);
        end
        
        
        if(r<=pk)
            a(s2+1,k+1) = -a(s1+1,k-1+1)/N(pk+1+1,r+1);
            d = d + a(s2+1,k+1)*N(r+1,pk+1);
        end
        ders(k+1,r+1) = d;
        j = s1; s1 = s2; s2=j;
    end
end


% Multiply through by the correct factors

r = p;
for k = 1 : 2
    for j = 0:p
        ders(k+1,j+1) = ders(k+1,j+1)*r;
    end
    r = r * (p-k);
end



end