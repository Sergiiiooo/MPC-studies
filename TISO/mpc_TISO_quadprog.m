% Controle Adaptativo Preditivo
% José Sergio Cruz Dantas Junior

% Simulaçao do Poster 'Observer Based Approach for the Economic Predictive Control of a TISO System'
% por linha de codigo usando a funçao 'quadprog'

clear; clc; close all; warning off;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');
%% Modelo CARIMA
Ts = 24;   % tempo de amostragem

% G(z^-1) = z^-d * b(z^-1)/a(z^-1)

d1 = 0; % atrasos
d2 = 0; % <-
Gz1  = tf([0 0.1],[1 -0.99],Ts,'Variable','z^-1','IODelay',d1)  % modelo da temperatura
Gz2  = tf([0 -0.5],[1 -0.47],Ts,'Variable','z^-1','IODelay',d2) % modelo da umidade
[b1,a1] = tfdata(Gz1,'v')
[b2,a2] = tfdata(Gz2,'v')

% b(z^-1) = 0 + b1 z^-1 + b2 z^-2 + ...
% a(z^-1) = 1 + a1 z^-1 + a2 z^-2 + ...
% c(z^-1) = 1 + c1 z^-1 + c2 z^-2 + ...

c = [1 -0]      % [1 c1 c2 ...]
Delta = tf([1 -1],[1],Ts,'Variable','z^-1') % Δ = 1 - z^-1
[delta, ~] = tfdata(Delta, 'v');

% a_til*y(k) = b1_til*Δu1(k) + b2_til*u2(k) + c*e(k)

a_til = conv(conv(delta,a1),a2) % a_til = Δ*a1*a2
b1_til = [zeros(1,d1) conv(a2,b1)]  % b1_til = a2*b1
b2_til = [zeros(1,d2) conv(conv(delta,a1),b2)]  % b2_til = Δ*a1*b2

% ajustando os vetores para terem o mesmo tamanho
Max = max([length(a_til)-1,length(b1_til)-1,length(b2_til)-1,length(c)-1])
a_til = [a_til zeros(1,Max-(length(a_til)-1))]
b1_til = [b1_til zeros(1,Max-(length(b1_til)-1))]
b2_til = [b2_til zeros(1,Max-(length(b2_til)-1))]
c = [c zeros(1,Max-(length(c)-1))]

% return
%% Espaço de Estados na forma canônica observável

% x(k+1) = A*x(k) + B1*Δu1(k) + B2*u2(k) + D*e(k)
% y(k) = H*x(k) + e(k)

%     ___                               ___          ___          ___
%     | -a_til(1)     1     0   ...     0 |          |  bi_til(1)   |
%     | -a_til(2)     0     1   ...     0 |          |  bi_til(2)   |
% A = |                 ...               |     Bi = |    ...       |
%     | -a_til(n-1)   0     0   ...     1 |          |  bi_til(n-1) |
%     | -a_til(n)     0     0   ...     0 |          |  bi_til(n)   |
%     ___                               ___ nxn      ___          ___ nx1
%     ___               ___
%     | c(1)-a_til(1)     |
%     | c(2)-a_til(2)     |          ___                ___
% D = |      ...          |      H = | 1    0   ...     0 |
%     | c(n-1)-a_til(n-1) |          ___                ___ 1xn
%     | c(n)-a_til(n)     |
%     ___               ___ nx1

A = [-a_til(2:end)' [eye(length(a_til)-2); zeros(1,length(a_til)-2)]]

B1 = [b1_til(2:end)']

B2 = [b2_til(2:end)']

D = [c(2:end)'-a_til(2:end)']

H = [1 zeros(1,Max-1)]

% return
%% Predições
% Horizonte de predição
N = 4

% Horizonte de controle
Nu = 2

% y = G1*Δu1 + G2*u2 + F*x + E*e

%      ___                                        ___
%      | H*Bi           0               ...      0  |
%      | H*A*Bi         H*Bi            ...      0  |
% Gi = | H*A^2*Bi       H*A*Bi          ...      0  |
%      |                 ...                        |
%      | H*A^(N-1)*B    H*A^(N-2)*B     ...     H*B |
%      ___                                        ___ NxNu
for i=1:N
   G1(i,:) = H*A^(i-1)*B1;
end
for j=2:Nu
    G1 = [G1 [0; G1(1:end-1,j-1)]];
end

%      ___                        ___            ___               ___
%      | h1    0     0     0    ... |    Nu=2    | h1   0            |
%      | h2    h1    0     0    ... |       \    | h2   h1           |
% G2 = | h3    h2    h1    0    ... |    ----|   | h3   h2+h1        |      
%      | h4    h3    h2    h1   ... |       /    | h4   h3+h2+h1     |
%      |          ...               |            |       ...         |
%      ___                        ___ NxN        ___               ___ NxNu
for i=1:N
   G2(i,:) = H*A^(i-1)*B2;
end
for j=2:N
    G2 = [G2 [0; G2(1:end-1,j-1)]];
end
for k=Nu+1:N
    G2(:,Nu) = G2(:,Nu) + G2(:,k);
end
G2 = resize(G2,[N,Nu]);

%      ___   ___
%      |  H*A  |
%      | H*A^2 |
%  F = | H*A^3 |
%      |  ...  |
%      | H*A^N |
%      ___   ___ Nxn
for i=1:N
    F(i,:) = H*A^i;
end
%      ___         ___
%      |     H*D     |
%      |    H*A*D    |
%  E = |   H*A^2*D   |
%      |    ...      |
%      | H*A^(N-1)*D |
%      ___         ___ Nx1
for i=1:N
    E(i,1) = H*A^(i-1)*D;
end

% G
G = [G1 G2]
F
E

% return
%% Simulacao
% parametros
qntd_k = 75*3;    % quantidade de amostras      75*24s = 1800s = 30min
tsim = qntd_k*Ts;   % tempo de simulaçao

cond_inicial = 0;
cond_desej = 1;

r = cond_desej - cond_inicial;  % referencia
t_r = 10*Ts;    % tempo de aplicaçao do sinal de referencia
ur = 1;         % segunda entrada
t_ur = ceil(qntd_k/2)*Ts;   % tempo de aplicaçao do sinal de ur

ref = r*ones(qntd_k,1); % vetor de referencia
Ur = [zeros(ceil(t_ur/Ts),1); ur*ones(ceil(tsim/Ts)-ceil(t_ur/Ts),1)];  % vetor da segunda entrada
% q = [zeros(ceil(Tdist/Ts),1); DistAmplitude*ones(ceil((tsim-Tdist)/Ts),1)];

% Matrizes de ponderaçao
Q = 1 * eye(N,N);
R1 = 1 * eye(Nu,Nu);
% R1(1,1) = 0;
R2 = 1 * eye(Nu,Nu);
R = [R1 zeros(Nu); zeros(Nu) R2];

% condiçoes iniciais
f = zeros(N,1);
W = zeros(N,1);
y1 = 0;
y2 = 0;
du0 = 0;
u0 = 0;
u_1 = 0;
x0 = zeros(length(A), 1);
e0 = 0;
[~,z1] = filter(b1,a1,0);
[~,z2] = filter(b2,a2,0);

for k=1:ceil(tsim/Ts)

    y(k) = y1+y2;
    x(:,k) = A*x0 + B1*du0 + B2*u0 + D*e0;  % x(k) = Ax(k-1) + B1Δu1(k-1) + B2u2(k-1) + De(k-1)
    e(k) = y(k) - H*x(:,k); % e(k)= y(k) - \hat{y}(k)
    
    % resposata livre
    for i=1:N
       f(i,1) = F(i,:)*x(:,k) + E(i,1)*e(k); 
    end
    
    % ref
    for i=1:N
        W(i,1) = r; % case the reference is known
    end
    
    % equacao quadrtica HH
    HH = 2*(G'*Q*G + R);
    
    % equacao quadrtica b
    b = 2*(G'*Q*(f-W) - [zeros(Nu,1); R2*(ones(Nu,1)*Ur(k))]);

    % constraints
    du1_max = 0.3*ones(Nu,1);
    du1_min = -0.1*ones(Nu,1);
    u2_max = 1.2*ones(Nu,1);
    u2_min = 0*ones(Nu,1);

    % constraints' matrices     A.U <= B
    Acon = [eye(2*Nu,2*Nu); -eye(2*Nu,2*Nu)];
    Bcon = [du1_max; u2_max; -du1_min; -u2_min];
    
    % resolvedor
    opt = optimoptions('quadprog','Display', 'off');
    % Acon = []; Bcon = [];
    sol = quadprog(HH,b,Acon,Bcon,[],[],[],[],[],opt);  % J = Δu'HΔu + b'u  -> programação quadrática

    % sol -> U = [Δu1(1)   Δu1(2)   ...   Δu1(Nu)   u2(1)   u2(2)   ...   u2(Nu)]'

    du1(k) = sol(1);
    u1(k) = u_1 + du1(k);
    u2(k) = sol(Nu+1);

    [y1,z1] = filter(b1(2:end),a1,u1(k),z1);
    [y2,z2] = filter(b2(2:end),a2,u2(k),z2);

    du0 = du1(k);
    u_1 = u1(k);
    u0 = u2(k);
    x0 = x(:,k);
    e0 = e(k);
end

t = 0:Ts:tsim-Ts;
t = t/60;   % tempo em minutos
ref = ref+cond_inicial;   % condiçao inicial
y = y+cond_inicial;   % condiçao inicial

figure
subplot(3,1,1)
hold on
stairs(t,y, 'b', 'linewidth', 2)
stairs(t,ref, '--k', 'linewidth', 2)
grid on
% axis tight
xlabel('Time (min)')
ylabel('Output')
legend('y(t)','r(t)', 'location', 'best')
title('Output Signal')
subplot(3,1,2)
hold on
stairs(t,u1, 'b', 'linewidth', 2)
grid on
% axis tight
xlabel('Time (min)')
ylabel('$u_1$')
title('Temperature Control Signal')
subplot(3,1,3)
hold on
stairs(t,u2, 'r', 'linewidth', 2)
grid on
% axis tight
xlabel('Time (min)')
ylabel('$u_2$')
title('Humidity Control Signal')
sgtitle("\textbf{Output and Control Signals for N = }" + N + "\textbf{ and Nu = }" + Nu)

% fim