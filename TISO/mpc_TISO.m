% Controle Adaptativo Preditivo
% José Sergio Cruz Dantas Junior

% Simulaçao do Poster 'Observer Based Approach for the Economic Predictive Control of a TISO System'

clear; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex'); % adicionando a biblioteca do latex ao matlab
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
Max = max([length(a_til)-1,length(b1_til)-1,length(b2_til)-1,length(c)-1]);
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
%% Matrizes de Prediçao

% Y = G1 Δu1 + G2 u2 + f    ;    f = Fx + Ee

% Horizonte de predição
N = 4

% Horizonte de controle
Nu = 2

%      ___                                         ___
%      | H*B1           0               ...      0   |
%      | H*A*B1         H*B1            ...      0   |
% G1 = | H*A^2*B1       H*A*B1          ...      0   |
%      |                 ...                         |
%      | H*A^(N-1)*B1   H*A^(N-2)*B1    ...     H*B1 |
%      ___                                         ___ NxNu
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

G = [G1 G2]
F
E

% return
%% Controle
% Matrizes de ponderaçao
Q = 1 * eye(N,N);
R1 = 1 * eye(Nu,Nu);
% R1(1,1) = 0;
R2 = 1 * eye(Nu,Nu);
R = [R1 zeros(Nu); zeros(Nu) R2];

% Sintonia dos controladores        U = - H^(-1) b

% K1 => primeira linha de (H^-1)*G'*Q
K1 = inv(G'*Q*G + R)*G'*Q;
K1 = K1(1,:);

% K2 => linha Nu+1 de (H^-1)*G'*Q
K2 = inv(G'*Q*G + R)*G'*Q;
K2 = K2(Nu+1,:);

K1r = K1*ones(N,1);
K2r = K2*ones(N,1);

% K3_ => primeira linha de (H^-1)*Kbr
% K3 => linha Nu+1 de (H^-1)*Kbr
% Kbr => [zeros(Nu,1); R2*ones(Nu,1)*ur]
K3 = inv(G'*Q*G + R);
K3_ = K3(1,:);
K3 = K3(Nu+1,:);
K3_ = K3_*[zeros(Nu,1); R2*ones(Nu,1)*1];
K3 = K3*[zeros(Nu,1); R2*ones(Nu,1)*1];

% Matrizes de controle
Kr = [K1r; K2r]
K = [K1*F; K2*F]
Ke = [K1*E; K2*E]
Kur = [K3_; K3]

% return
%% Simulação
qntd_k = 75*3;    % quantidade de amostras      75*24s = 1800s = 30min
tsim = qntd_k*Ts;   % tempo de simulaçao

cond_inicial = 0;
cond_desej = 1;

r = cond_desej - cond_inicial;  % referencia
t_r = 0*Ts;    % tempo de aplicaçao do sinal de referencia
ur = 1;         % segunda entrada
t_ur = ceil(qntd_k/2)*Ts;   % tempo de aplicaçao do sinal de u2r

out = sim('mpc_2input_simulink.slx');

t = out.t/60;   % tempo de plotagem em minutos
ref = out.r+cond_inicial;   % condiçao inicial
y = out.y+cond_inicial;   % condiçao inicial

figure
subplot(3,1,1)
hold on
stairs(t,y,'b','LineWidth',2)
stairs(t,ref, 'k--','LineWidth',2)
% plot(t,out.ur, 'r--')
grid on
% axis tight
xlabel('Time (min)')
ylabel('Output')
% legend('y(t)','r(t)','ur(t)', 'location', 'best')
legend('y(t)','r(t)', 'location', 'best')
title('Output Signal')
subplot(3,1,2)
hold on
stairs(t,out.u1,'b','LineWidth',2)
grid on
% axis tight
xlabel('Time (min)')
ylabel('$u_1$')
title('Temperature Control Signal')
subplot(3,1,3)
hold on
stairs(t,out.u2,'r','LineWidth',2)
grid on
% axis tight
xlabel('Time (min)')
ylabel('$u_2$')
title('Humidity Control Signal')
sgtitle("\textbf{Output and Control Signals for N = }" + N + "\textbf{ and Nu = }" + Nu)

% fim