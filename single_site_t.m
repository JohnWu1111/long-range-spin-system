clear;
% close all;
clc;
format long
tic;

% Definition of parameters
N = 2000; %size
lN = length(N);
J = 1;
h = 0.3;
T = 0:1:1000;

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
S1_y = (S1_p - S1_m)/2i;

% construction of Hamiltonian
H1 = -J*S1_z.^2/(2*N1);
H2 = h*S1_x;
H = diag(H1 + J/8) + H2;

% time revolution
[V,D] = eig(H);
e = diag(D);

phi0 = zeros(N+1,1);
phi0(1) = 1;

temp = V'*phi0;
trans = exp(-1i*e*T);
temp = trans.*temp;
phit = V*temp;
phi_norm = sum(abs(phit).^2);

% mz = sum(conj(phit).*(S1_z.*phit))/N;
mz = sum(S1_z.*abs(phit).^2)/N;
mx = real(sum(conj(phit).*(S1_x*phit)))/N;
my = real(sum(conj(phit).*(S1_y*phit)))/N;
m = sqrt(mx.^2+my.^2+mz.^2);
mz0 = sum(S1_z.*abs(V(:,1)).^2)/N;

figure;
plot(T,mz)

len = length(e);
e = sort(e);
s = zeros(len-1,1);
r = zeros(len-2,1);
for i = 1:len-1
    s(i) = e(i+1) - e(i);
end

for i = 1:len-2
    r(i) = s(i+1)/s(i);
end

q = min([r 1./r],[],2);

% figure;
% histogram(q,100,'Normalization','pdf','DisplayStyle','stairs');
% 
% mean(q)

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