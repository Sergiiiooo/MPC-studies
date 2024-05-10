% Controle Adaptativo Preditivo
% Exemplo 2 - Caso 2
% José Sergio Cruz Dantas Junior

clc; clear; close all; warning off;

%% Dados iniciais

b = [0 -0.2 1];
a = [1 -1 0];
c = [1 -2.7 2.43 -0.729]; %(1 - 0.9)^3

% Horizonte de predição
N = 70

% Horizonte de Controle
Nu = 1

% Q
Q = 1*eye(N)

% R
R = 0*eye(Nu)

%% Transfer Function

numGs = b;
denGs = a;
Ts = 0.01;
% s = tf('s');
z = tf('z', Ts);

display("Função de Transferência:")
Gs = tf([numGs],[denGs],'Variable','s','IODelay',0)
Gz = zpk(c2d(Gs,0.01,'zoh'))

[b,a] = tfdata(Gz,'v');

% integral
Delta = 1/(1 - z^-1);
[numDelta, denDelta] = tfdata(Delta, 'v');

% a ~
a_ = conv(a, denDelta);

% correção do c
c_ = zeros(1, length(a_));
for i=length(c)+1:length(a_)
    c(i) = 0;
end

%% State-Space model

Max = max([length(b),length(a_)-1,length(c)-1]);

for i=1:Max
    for j=1:Max
        if j-1==i
            A(i,j) = 1;
        else if j-1~=i
                A(i,j) = 0;
        end
        end
        if j==1
            if i>length(a_)-1
                A(i,j) = 0;
            else
                A(i,j) = -a_(1,i+1);
            end
        end
    end
end

for i=1:Max
    if i>length(b)
                B(i,1) = 0;
            else
                B(i,1) = b(1,i);
            end
end

for i=1:Max
    if i>length(c)-1
        if i>length(a_)-1
            D(i,1) = 0;
        else
            D(i,1) = -a_(1,i+1);
        end
    else
        if i>length(a_)-1
            D(i,1) = c(1,i+1);
        else
            D(i,1) = c(1,i+1)-a_(1,i+1);
        end
    end
end

for i=1:Max
    if i==1
        H(1,i) = 1;
    else
        H(1,i) = 0;
    end
end

display("Modelo de Espaço de Estados:")
A
B
D
H

%% Matrizes de Predição

% Matriz G
G = zeros(N,Nu);

for i=1:N
   G(i,1) = H*A^(i-1)*B;
end

for i=2:N
   for j=2:Nu
      G(i,j) = G(i-1,j-1);
   end
end

% Matriz F
F = zeros(N,length(a_)-1);

for i=1:N
   F(i,:) = H*A^i; 
end

% Matriz E
E = zeros(N,1);

for i=1:N
   E(i,1) = H*A^(i-1)*D; 
end

% K
K = inv(G'*Q*G + R)*G'*Q;
K = K(1,:);

% Kr
Kr = sum(K);

%% Simulação

Gsr = Gz;

% Parâmetros
Tsim = 20;
Tref = 0;
Tdist = 10;
RefAmplitude = 1;
RefInitialValue = 0;
DistAmplitude = -1;

sim = sim('simulink_MPC2');

figure(1)
subplot(2,1,1)
plot(sim.t,sim.y,'b', linewidth=1.5)
axis ([0 20 -0.2 1.2])
grid on
ylabel('y(t)')

subplot(2,1,2)
plot(sim.t,sim.u,'b', linewidth=1.5)
axis ([0 20 -2 6])
grid on
ylabel('u(t)')
xlabel('time(s)')
