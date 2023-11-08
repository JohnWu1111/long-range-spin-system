clear;
% close all;
clc;
format long
tic;

myseed = 1;
rng(myseed);

%% Definition of parameters
L = 100; %size
J = 1;
h = 0.4;
dt = 1e-3;
t = 0:dt:1000;
nt = length(t);
alpha = 100;
T = 0.1;
lambda = 1;
D = sqrt(2*lambda*T);
x = 1:L;

S = zeros(3,L);
% S(1,:) = -0.11;
% S(3,:) = sqrt(1/4-S(1,:).^2);
% mz = S(3,:);
S(3,:) = 1;
mz =1;
mz_t = zeros(1,nt);
mz_t(1) = mz;
L_S = zeros(L,nt);
L_S(:,1) = 1;

dist = min(abs(x - x'),abs(x - (x+L)'));
dist = min(dist,abs(x+L - x'));
coeff = 1./(dist.^alpha);
coeff(coeff==Inf) = 0;
xx = 1:L/2;
fL = sum(1./xx.^alpha);

%% time evolution
for i = 2:nt
    dS = zeros(3,L);
    noise = D*randn(3,L)/sqrt(dt);
    H = noise;
    H(1,:) = H(1,:) + h;
%     H(3,:) = H(3,:)-J*sum(S(3,:).*coeff)/fL;
    H(3,:) = H(3,:)-J*S(3,:)*coeff/fL;
    HcX = cross(H,S);
    XcHcX = cross(S,HcX);
    S_temp = S + (HcX - lambda*XcHcX)*dt;
    H_temp = noise;
    H_temp(1,:) = H_temp(1,:) + h;
    H_temp(3,:) = H_temp(3,:)-J*S_temp(3,:)*coeff/fL;
    HcX_temp = cross(H_temp,S_temp);
    XcHcX_temp = cross(S_temp,HcX_temp);
    S = S + (HcX - lambda*XcHcX + HcX_temp - lambda*XcHcX_temp)*dt/2;
    L_S(:,i) = sqrt(sum(S.^2))';
    S = S./L_S(:,i)';
    mz = sum(S(3,:))/L;
    mz_t(i) = mz;    
end

filename = strcat('L = ',num2str(L), ', h = ', num2str(h),', lambda = ', num2str(lambda), ', T = ', num2str(T), ', alpha = ', num2str(alpha), ', myseed = ', num2str(myseed));
figure('Name',filename);
set(gcf, 'position', [250 70 1400 900]);
plot(t,mz_t)

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