%% Activity nº 4
% Model Predictive Control
clc; clear; close all; warning off;

% Parameters
b = [0 0.2713]
a = [1 -0.8351]
c = [1 0 0]

N = 10  % Prediction Horizon
Nu = N  % Control Horizon
d = 4;  % I/O Delay

% Ponderation matrices
Q = 1*eye(N);
R = 100*eye(Nu);

Ts = 0.1    % Sample Time
Gz = tf([b],[a],Ts,'Variable','z^-1','IODelay',d)   % Process
Delta = tf([1],[1 -1],Ts,'Variable','z^-1','IODelay',0) % Delta

% a ~
a_til = conv(a, [1 -1])

% adjusting the vectors to have the same length
Max = max([length(a_til)-1,length(b)-1,length(c)-1]);
a_til = [a_til zeros(1,Max-(length(a_til)-1))];
a = [a zeros(1,Max-(length(a)-1))];
b = [b zeros(1,Max-(length(b)-1))];
c = [c zeros(1,Max-(length(c)-1))];

%% State-Space formulation observable cannonical form
A = [-a_til(2:end)' [eye(length(a_til)-2); zeros(1,length(a_til)-2)]]
B = [b(2:end)']
D = [c(2:end)'-a_til(2:end)']
H = [1 zeros(1,Max-1)]

%% Prediction Matrices
for i=1:N
   G(i,1) = H*A^(i-1)*B;
end
G = toeplitz(G);
G = tril(G)

for i=1:N
   F(i,:) = H*A^i; 
end
F

for i=1:N
   E(i,1) = H*A^(i-1)*D; 
end
%% Control
K = inv(G'*Q*G + R)*G'*Q;
K = K(1,:)
Kr = sum(K)
KF = K*F

%% Simulation
Tsim = 30;  % Simulation time
Tdist = Tsim/2; % Disturbance time
ref = 1;    % Reference Amplitude
dist = -0.5;   % Disturbance Amplitude

sim = sim('simu_MPC_deadtime');

figure
subplot(3,1,1)
stairs(sim.y,'b', linewidth=1.5)
xlabel('discrete time')
ylabel('y(k)')
axis("tight")
grid on
subplot(3,1,2)
stairs(sim.u,'r', linewidth=1)
xlabel('discrete time')
ylabel('u(k)')
axis("tight")
grid on
subplot(3,1,3)
stairs(sim.du,'r', linewidth=1)
xlabel('discrete time')
ylabel('Δu(k)')
axis("tight")
grid on

fprintf('\n\n Code Author: Jose Sergio Cruz Dantas Junior - UFC \n')
% fim