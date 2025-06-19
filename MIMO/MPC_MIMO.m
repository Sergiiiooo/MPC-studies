%% Model Predictive Control for MIMO systems
% Example of the page 139 of Camacho's Book
clc; clear; close all; warning off;

% Process
G1 = tf(1,[0.7 1]);
G2 = tf(5,[0.3 1]);
G3 = tf(1,[0.5 1]);
G4 = tf(2,[0.4 1]);
GG = [G1 G2; G3 G4]

Ts = 0.03;  % Sample Time
GGd = c2d(GG,Ts,'zoh')  % Discretization

[b1,a1] = tfdata(GGd(1,1),'v');
[b2,a2] = tfdata(GGd(1,2),'v');
[b3,a3] = tfdata(GGd(2,1),'v');
[b4,a4] = tfdata(GGd(2,2),'v');
c1 = [1 0 0];
c2 = [1 0 0];

% a ~
delta = [1 -1];
Delta = tf([1],delta,Ts,'Variable','z^-1','IODelay',0); % Delta
Delta = [Delta 0; 0 Delta];
a1t = conv(delta,conv(a1,a2));
a2t = conv(delta,conv(a3,a4));
b1t = conv(b1,a2);
b2t = conv(b2,a1);
b3t = conv(b3,a4);
b4t = conv(b4,a3);

% adjusting the vectors to have the same length
Max = max([length(a1t)-1,length(a2t)-1,length(b1t)-1,length(b2t)-1,length(b3t)-1,length(b4t)-1,length(c1)-1,length(c2)-1]);
a1t = [a1t zeros(1,Max-(length(a1t)-1))];
a2t = [a2t zeros(1,Max-(length(a2t)-1))];
b1t = [b1t zeros(1,Max-(length(b1t)-1))];
b2t = [b2t zeros(1,Max-(length(b2t)-1))];
b3t = [b3t zeros(1,Max-(length(b3t)-1))];
b4t = [b4t zeros(1,Max-(length(b4t)-1))];
c1 = [c1 zeros(1,Max-(length(c1)-1))];
c2 = [c2 zeros(1,Max-(length(c2)-1))];

N = 3  % Prediction Horizon
Nu = 2  % Control Horizon

% Ponderation matrices
Q = 1*eye(2*N);
R = 0.05*eye(2*Nu);

%% State-Space formulation observable cannonical form
% Augmented Matrices
A1 = [-a1t(2:end)' [eye(length(a1t)-2); zeros(1,length(a1t)-2)]];
A2 = [-a2t(2:end)' [eye(length(a2t)-2); zeros(1,length(a2t)-2)]];
A = [A1 zeros(size(A1)); zeros(size(A2)) A2]
B1=[b1t(2:end)]';
B2=[b2t(2:end)]';
B3=[b3t(2:end)]';
B4=[b4t(2:end)]';
B=[B1 B2; B3 B4];
D1 = [c1(2:end)'-a1t(2:end)'];
D2 = [c2(2:end)'-a2t(2:end)'];
D = [D1 zeros(size(D1)); zeros(size(D2)) D2]
H = [1 zeros(1,Max-1)];
H = [H zeros(size(H)); zeros(size(H)) H]

%% Prediction Matrices
G = zeros(2*N,2*Nu);
for i=1:N
    for j=1:Nu
        if i>=j
            G(2*i-1:2*i,2*j-1:2*j) = H*A^(i-j)*B;
        end
    end
end
G

F = zeros(N * size(H,1), size(A,2));
for i = 1:N
    row_idx = (i-1)*size(H,1)+1 : i*size(H,1);
    F(row_idx, :) = H * A^i;
end
F

E = zeros(N * size(H,1), size(D,2));
for i = 1:N
    row_idx = (i-1)*size(H,1)+1 : i*size(H,1);
    E(row_idx, :) = H * A^(i-1) * D;
end
E

%% Control
K = inv(G'*Q*G + R)*G'*Q;
K = K(1:2,:)
KF = K*F
KE = K*E
Kr = [sum(K(1, 1:2:2*N)), sum(K(1, 2:2:2*N));
      sum(K(2, 1:2:2*N)), sum(K(2, 2:2:2*N))];

%% Simulation
Tsim = 100*Ts;
r1 = 0.5;
r1d = -0.1;
Tr1d = 50*Ts;
r2 = 0.3;

sim = sim('simu_MPC_MIMO');

figure(1)
plot(sim.y(:,1),'r', linewidth=1.5)
hold on
plot(sim.y(:,2),'b', linewidth=1.5)
plot(sim.ref(:,1),'--k', linewidth=1)
plot(sim.ref(:,2),'--k', linewidth=1)
grid on
axis ([0 100 0 0.6])
xlabel('Samples')
ylabel('Outputs and References')
legend('y1', 'y2', 'location', 'best')

figure(2)
plot(sim.u(:,1),'r', linewidth=1.5)
hold on
plot(sim.u(:,2),'b', linewidth=1.5)
grid on
axis ([0 100 -0.2 0.8])
xlabel('Samples')
ylabel('Inputs')
legend('u1', 'u2', 'location', 'best')

fprintf('\n\n Code Author: Jose Sergio Cruz Dantas Junior - UFC \n')