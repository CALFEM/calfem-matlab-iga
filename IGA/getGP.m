function [ gp,w ] = getGP( p )
%[ gp,w ] = getGP( p )
%-------------------------------------------------------------
% PURPOSE:
%   Given the order p returns appropriate Gauss-points
%   and weghts for Gauss-Legandre quadrature.
%
% INPUT: p = polynomial order
%
% OUTPUT: [gp,w] = univariate Gauss-points and weights for
%                  parent domain in [-1,1]
%-------------------------------------------------------------

n = ceil((p*2+1)/2);

if n == 2
    gp = [-1/sqrt(3) 1/sqrt(3)];   % Two gp's
    w = [1 1];                     %
elseif n == 3
    w = [5/9 8/9 5/9];                 % Three gp's
    gp = [-sqrt(3/5) .0 sqrt(3/5)];    %
elseif n == 4
    w = [.347854845137454 .652145154862546 .652145154862546 .347854845137454];      % Four gp's
    gp = [-.861136311594053 -.339981043584856 .339981043584856 .861136311594053];   %
elseif n == 5
    w = [(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];
    gp = [-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))];
elseif n == 6
    w = [.171324492379170 .360761573048139 .467913934572691 .467913934572691 .360761573048139 .171324492379170];
    gp = [-.932469514203152 -.661209386466265 -.238619186083197 .238619186083197 .661209386466265 .932469514203152];
end
end

