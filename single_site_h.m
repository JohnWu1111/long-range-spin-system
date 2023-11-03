clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 5000; %size
lN = length(N);
J = 1;
h = [0.1 0.2 0.4 0.6 0.8 1 1.2 1.5 1.8 2];
lh = length(h);

N1 = N;
S1 = N1/2;

S1_z = zeros(N1+1,1);
S1_p = zeros(N1+1);
S1_m = zeros(N1+1);

% construction of matrice
for m = 1:N1+1
    S1_z(m) = N1/2 - (m-1);
end

for m = 1:N1
    S1_p(m,m+1) = sqrt(S1*(S1+1)-S1_z(m+1)*(S1_z(m+1)+1));
    S1_m(m+1,m) = sqrt(S1*(S1+1)-S1_z(m)*(S1_z(m)-1));
end

S1_x = (S1_p + S1_m)/2;

mz = zeros(1,lh);
mx = zeros(1,lh);

for n = 1:lh
    
    % construction of Hamiltonian
    H1 = -J*S1_z.^2/(2*N1);
    H1 = diag(H1);
    H2 = h(n)*S1_x;
    H = H1 + H2;
    
    H = sparse(H);
    
    % time revolution
    [V,D] = eigs(H,2,'smallestreal');
    e = diag(D);
    V = V(:,1);
    
    Mz2 = S1_z.^2/N;
    Mx = S1_x/N;
    mz(n) = Mz2'*(V.^2);
    mx(n) = V'*Mx*V;
    
end

figure;
plot(h,mz)
xlabel('h')
ylabel('S_z^2/N')

figure;
plot(h,mx)
xlabel('h')
ylabel('S_x/N')

toc;

function y = kron4(a,b,c,d)
y = kron(kron(kron(a,b),c),d);
end

function y = kron3(a,b,c)
y = (kron(kron(a,b),c));
end

function y = kron_p(a,b)
la = length(a);
lb = length(b);
y = zeros(la*lb,1);
for i = 1:la
    for j = 1:lb
        y((i-1)*lb+j) = a(i) + b(j);
    end
end
end

function y = kron_p4(a,b,c,d)
y = kron_p(kron_p(kron_p(a,b),c),d);
end