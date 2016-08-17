function [Xibar,Qw] = RefineKnotVectCurve(n,p,Xi,Pw,X,r)
%--------------------------------------------------------------
%function [Ubar,Qw] = RefineKnotVectCurve(n,p,U,Pw,X,r)
% algorithm A5.4 in The NURBS Book (Piegl & Tiller)
% insert multiple knots into curve
%INPUT:
% n         : number ob basis functions -1 !
%        NURBS-Book: n+1 # basis, np max index (startindex 0)
%        here        n   # basis and max index (startindex 1)
% p          : degree of the basis functions
% Xi         : old knotvector
% Pw         : old control points
% X          : vector of new knots (multiple entries possible)
% r          :  (size of X) -1 (count the multple entries as well
%                reason: same as above: max X index
%OUTPUT:
% Xibar       : new knot vector
% Qw         : new control points
%--------------------------------------------------------------

%initialise arrays;
dim  = size(Pw,2);
Qw   = zeros(n+r+2,dim);
Xibar = zeros(1,n+p+1+r);
%
m = n+p+1;
a = FindSpan(n,p,X(1),Xi);
b = FindSpan(n,p,X(r+1),Xi);
b = b+1;

for j=0:a-p
    Qw(j+1,:) = Pw(j+1,:);
end
for j=b-1:n
    Qw(j+r+2,:) = Pw(j+1,:);
end
for j=0:a
    Xibar(j+1)= Xi(j+1);
end
for j=b+p :m
    Xibar(j+r+2) = Xi(j+1);
end
i=b+p-1;
k=b+p+r;
for j=r:-1:0
    while (X(j+1) <= Xi(i+1) && i>a)
        Qw(k-p,:) = Pw(i-p,:);
        Xibar(k+1) = Xi(i+1);
        k=k-1;
        i=i-1;
    end
    Qw(k-p,:) = Qw(k-p+1,:);
    for l=1:p
        ind = k-p+l;
        alfa = Xibar(k+l+1) - X(j+1);
        if (abs(alfa) == 0)
            Qw(ind,:) = Qw(ind+1,:);
        else
            alfa = alfa/(Xibar(k+l+1) - Xi(i-p+l+1));
            Qw(ind,:) = alfa* Qw(ind,:) + (1-alfa)* Qw(ind+1,:);
        end
    end
    Xibar(k+1) = X(j+1);
    k=k-1;
end
end