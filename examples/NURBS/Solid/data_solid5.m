p = 2;
q = 1;
r = 2;

Xi = [0,0,.5,1,1];
Xi = [Xi(1) Xi Xi(end)];
Eta = [0,1];
Eta = [Eta(1) Eta Eta(end)];
Zeta = [0,1];
Zeta = [Zeta(1) Zeta Zeta(end)];
k_ = length(Xi);
l_ = length(Eta);
m_ = length(Zeta);


B = cell(4,2,2);

s = 1/sqrt(2);

B{1,1,1} = [0 -1 -1 1];
B{1,1,2} = [0 -1 1 1];
B{1,2,1} = [0 1 -1 1];
B{1,2,2} = [0 1 1 1];

B{2,1,1} = [2 -1 0 1];
B{2,1,2} = [2 -1 2 1];
B{2,2,1} = [2 1 0 1];
B{2,2,2} = [2 1 2 1];

B{3,1,1} = [4 -1 -1 1];
B{3,1,2} = [4 -1 1 1];
B{3,2,1} = [4 1 -1 1];
B{3,2,2} = [4 1 1 1];

B{4,1,1} = [6 -1 -1 1];
B{4,1,2} = [6 -1 1 1];
B{4,2,1} = [6 1 -1 1];
B{4,2,2} = [6 1 1 1];


th=pi/4;
R = [ 1 0 0 0; 0 cos(th) -sin(th) 0; 0 sin(th) cos(th) 0; 0 0 0 1];
th=pi/4/3;
R2 = [ cos(th) 0 -sin(th) 0; 
       0 1 0 0; 
       sin(th) 0 cos(th) 0; 
       0 0 0 1];
for j = 1 : size(B,2)
    for k = 1 : size(B,3)
        B{2,j,k} = (R2*(R*B{2,j,k}'))';
    end
end

for j = 1 : size(B,2)
    for k = 1 : size(B,3)
        B{3,j,k} = (R2*R2*R*R*B{3,j,k}')';
    end
end

for j = 1 : size(B,2)
    for k = 1 : size(B,3)
        B{4,j,k} = (R2*R2*R2*R*R*R*B{4,j,k}')';
    end
end

% Number of weights
n = size(B,1);
m = size(B,2);
l = size(B,3);

% Order of basis
p = k_-n-1;
q = l_-m-1;
r = m_ -l -1;
