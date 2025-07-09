% MPC por Eq. Diofantina
% Ex 3.3 - Camacho
clear; clc

A = [1 -0.8]
na = length(A);
delta = [1 -1];
Atil = conv(delta,A);
N = 5
B = [0.4 0.6]
lambda = 0.8;
R = lambda*eye(N);

P1 = [1 zeros(1,na)];
P2 = Atil;

F = [];
Gp = [];

for n=1:N
    [q,r] = deconv(P1,P2);
    F = [F; r(2:end)];
    P1 = [r(2:end) 0];
    E(n) = q;
    aux = conv(E,B);
    Gp = [Gp; aux(n+1:end)];
end

Ga = aux(1:n);

tic
for i=1:N
   G(i,:) = Ga(i);
end
for j=2:N
    G = [G [0; G(1:end-1,j-1)]];
end
Teste1 = toc
G
clear G
tic
G = zeros(N);
for i = 1:N
    G(i:size(G,1)+1:N^2) = Ga(i);
end
G = tril(G);
Teste2 = toc
G
clear G
tic
G = toeplitz(Ga);
G = tril(G);
Teste3 = toc
G

aux = inv(G'*G+R)*G';
K = aux(1,:);
Kr = sum(K)
KGp = [0 K*Gp]
KF = K*F