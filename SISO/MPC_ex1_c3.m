% Controle Adaptativo Preditivo
% Exemplo 1 - Caso 3
% José Sergio Cruz Dantas Junior

clc; clear; close all; warning off;

%% Dados iniciais

b = [0.4 0.6];
a = [1 -0.8];
c = [1 -1.6 0.64];

% Horizonte de predição
N = 3

% Horizonte de Controle
Nu = N

% Q
Q = 1*eye(N)

% R
R = 0.8*eye(Nu)

%% Transfer Function

numGz = b;
denGz = a;
Ts = 1;
% s = tf('s');
z = tf('z', Ts);

display("Função de Transferência:")
Gz = tf([numGz],[denGz],Ts,'Variable','q^-1','IODelay',1)

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

%% Indice de robustez

% Incerteza multiplicativa
wz = logspace(-1,log10(pi/Ts),1000);
dP = abs(freqresp((Gz - Gz*z^-1),wz))./abs(freqresp(Gz*z^-1,wz));
dP = squeeze(dP(1,1,:));

aux = log10(pi/Ts);
n = 1000;
w = logspace(-2,aux,n);
Ir = zeros(1,n);
for k = 1:n
    zw = freqresp(z,w(k));
Gzw = zw^(-1)*(0.4+0.6*zw^-1)/(1 - 0.8*zw^-1);
V = K*E + (K*F-K*E*H)*((zw*eye(size(A))-A+D*H )\D);
ir = (Gzw*V\((K*F)*((zw*eye(size(A))-A)\B)+1))*(1-zw^(-1));
Ir(:,k) = abs(ir);
end

%% Simulação

Gsr = Gz;

% Parâmetros
Tsim = 60;
Tref = 0;
Tdist = 20;
RefAmplitude = 1;
RefInitialValue = 0;
DistAmplitude = -0.1;

sim = sim('simulink_MPC');

figure(1)
plot(sim.y,'b', linewidth=1.5)
hold on
stairs(sim.du,'r', linewidth=1)
axis ([1 60 -0.4 1.2])
grid on
xlabel('discrete time')
legend('y(t)', 'Δu(t)', 'location', 'best')

figure(2)
loglog(w, dP, 'k', linewidth=1)
hold on
loglog(w, Ir, 'b', linewidth=1.5)
% grid on
axis ([1e-2 1e1 1e-2 1e3])
xlabel('frequency ω')
legend('Incerteza δP(ω)','Indice de robustez Ir(ω)','location', 'best')
