function [ INC,IEN ] = BldINCIEN( p,q,r,n,m,l )
%[ INN,IEN ] = BldINCIEN( p,q,r,n,m,l )
%[ INN,IEN ] = BldINCIEN( p,q,n,m )
%[ INN,IEN ] = BldINCIEN( p,n )
%-------------------------------------------------------------
% PURPOSE:
%  Constructs the INC and IEN arrays. 
%
% INPUT: p = polynomial order in first direction
%        q = polynomial order in second direction
%        q = polynomial order in third direction
%        n = number of control points in first direction
%        m = number of control points in second direction
%        l = number of control points in third direction
%
% OUTPUT: INC = NURBS coordinates (See Hughes et. al. 2009) (Relates a
%               global function number to uni-variate knot indices).
%               INC(A,:) = xi_coord(A),eta_coord(A),[zeta_coord(A)], where 
%               A is global basis function number, xi_coord(A) is the
%               xi-coordinate of A, etc.
%         IEN = Relates global basis function numbers to their local
%               orderings, often referred to as "element nodes" array.
%               IEN(e,a) = A, where e ís element number, a is local basis
%               function number, and A is global basis function number.
%-------------------------------------------------------------


% ---------------------- for Solids ------------------------- %
if nargin == 6
    nel = (n-p) * (m-q) * (l-r); % number of elements
    nnp = n*m*l; % number of global basis functions
    nen = (p+1)*(q+1)*(r+1); % number of local basis functions
    INC = zeros([nnp],[3]); % NURBS coordinates array
    IEN = zeros([nen],[nel]); % connectivity array ( = Enod' )
    
    % Local variable initializations:
    %e, A, B, b, i, j, k, iloc, jloc, kloc % should all be
    % initialized to zero
    A = 0; e=0;
    
    for k = 1 : l
        for j = 1 : m
            for i = 1 : n
                
                A = A + 1;
                INC(A,1) = i;
                INC(A,2) = j;
                INC(A,3) = k;
                
                if i >= (p+1) & j >= (q+1) & k >= (r+1)
                    e = e + 1;
                    for kloc = 0 : r
                        for jloc = 0 : q
                            for iloc = 0 : p
                                B_ = A - kloc*n*m - jloc*n - iloc; % global function number
                                b = kloc*(p+1)*(q+1) + jloc*(p+1)+ iloc + 1; % local function number
                                IEN(b,e) = B_; % assign connectivity
                            end
                        end
                    end
                end
            end
        end
    end
    
% ---------------------------- For surfaces --------------------------- %
elseif nargin == 4
    m = n; % Fix var names.
    n = r; %
    
    
    nel = (n-p) * (m-q); % number of elements
    nnp = n*m; % number of global basis functions
    nen = (p+1)*(q+1); % number of local basis functions
    INC = zeros(nnp,2); % NURBS coordinates array
    IEN = zeros(nen,nel); % connectivity array ( = Enod' )
    
    % Local variable initializations:
    A = 0; e=0;
    
    for j = 1 : m
        for i = 1 : n
            
            A = A + 1;
            INC(A,1) = i;
            INC(A,2) = j;
            
            if i >= (p+1) & j >= (q+1)
                e = e + 1;
                for jloc = 0 : q
                    for iloc = 0 : p
                        B_ = A - jloc*n - iloc; % global function number
                        b =  jloc*(p+1)+ iloc + 1; % local function number
                        IEN(b,e) = B_; % assign connectivity
                    end
                end
            end
        end
    end
    
    
% ---------------------------- For curves --------------------------- %
elseif nargin == 2
    n = q; % Fix var names.
    
    
    nel = (n-p); % number of elements
    nnp = n; % number of global basis functions
    nen = (p+1); % number of local basis functions
    INC = zeros(nnp,1); % NURBS coordinates array
    IEN = zeros(nen,nel); % connectivity array ( = Enod' )
    
    % Local variable initializations:
    A = 0; e=0;
    
    for i = 1 : n

        A = A + 1;
        INC(A,1) = i;

        if i >= (p+1)
            e = e + 1;
            for iloc = 0 : p
                B_ = A - iloc; % global function number
                b = iloc + 1; % local function number
                IEN(b,e) = B_; % assign connectivity
            end
        end
    end
    
% ----------------------- Else: ------------------------- %
else
    fprintf('Error, only 4 or 6 input arguments supported.\n')
    INC = nan;
    IEN = nan;
    return
end

end

