clear;
% close all;
clc;
format long
tic;

%% Definition of parameters
L = 1; %size
J = 1;
h = 0.1;
dt = 1e-3;
T = 0:dt:1000;
nT = length(T);

% S = zeros(3,L);
% S(3,:) = 1/2;
% mz = 1/2;
% mz_t = zeros(1,nT);
% mz_t(1) = mz;
% L_S = zeros(L,nT);
% L_S(:,1) = 1/2;
% for i = 2:nT
%     dS = zeros(3,L);
%     dS(1,:) = -J*mz*S(2,:)/2;
%     dS(2,:) = h*S(3,:)-J*mz*S(1,:)/2;
%     dS(3,:) = -h*S(2,:);
%     S = S + dS*dt;
%     L_S(:,i) = sqrt(sum(S.^2))';
%     S = S./(2*L_S(:,i))';
%     mz = sum(S(3,:))/L;
%     mz_t(i) = mz;    
% end

S = zeros(3,1);
% S(1,:) = -0.11;
% S(3,:) = sqrt(1/4-S(1,:).^2);
% mz = S(3,:);
S(3,:) = 1/2;
mz = 1/2;
mz_t = zeros(1,nT);
mz_t(1) = mz;
L_S = zeros(L,nT);
L_S(:,1) = 1/2;
for i = 2:nT
    dS = zeros(3,L);
    dS(1,:) = -J*mz*S(2,:);
    dS(2,:) = h*S(3,:)+J*mz*S(1,:);
    dS(3,:) = -h*S(2,:);
    S = S + dS*dt;
    L_S(:,i) = sqrt(sum(S.^2))';
    S = S./(2*L_S(:,i))';
    mz = sum(S(3,:))/L;
    mz_t(i) = mz;    
end

figure;
plot(T,mz_t)

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