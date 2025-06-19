%% Model Predictive Control for MIMO systems
% Example of the page 139 of Camacho's Book using Quadprog function

clear all; clc; close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');
%% Model
% Process
G1 = tf(1,[0.7 1]);
G2 = tf(5,[0.3 1]);
G3 = tf(1,[0.5 1]);
G4 = tf(2,[0.4 1]);
GG = [G1 G2; G3 G4]

% Sampling time
Ts = 0.03;

GGd = c2d(GG,Ts,'zoh')  % Discretization

% return
%% Controller tuning parameters
% Prediction horizon
N = 3;
% Control horizon
Nu = 2;

% Weight values
R = 0.05*eye(2*Nu);
Q = eye(2*N);

% Constraints
du_min = -inf;
du_max = inf;
u_min = -inf;
u_max = inf;
y_min = -inf;
y_max = inf;
% return
%% CARIMA model
delta = [1 -1];
Delta = tf(1,delta,Ts,'Variable','z^-1'); % Δ = 1 - z^-1

% Polynomials b(z^-1) and a(z^-1)
[b1,a1] = tfdata(GGd(1,1),'v');
[b2,a2] = tfdata(GGd(1,2),'v');
[b3,a3] = tfdata(GGd(2,1),'v');
[b4,a4] = tfdata(GGd(2,2),'v');

% Disturbance rejection polynomial c(z^-1) 
c1  = conv([1 0],[1 0]);
c2  = conv([1 0],[1 0]);

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
% return
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
% return
%% Prediction Matrix
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
% return
%% Simulation parameters
% Time parameters
Tsim = 100*Ts;
qntd_k = ceil(Tsim/Ts);

% Reference parameters
r1 = 0.5;
r2 = 0.3;

% % Reference and disturbance arrays
% ref = [zeros(Tr1,1); r1*ones(Tr2-Tr1,1); (r1+r2)*ones(qntd_k-Tr2,1)];
% q = [zeros(Tpert,1); pert*ones(qntd_k-Tpert,1)];

%% GPC
% Quadratic equation HH
HH = 2*(G'*Q*G + R);

% Constraints matrix Ac
% Acon = [G -ones(size(G));
%     -G -ones(size(G));
%     ones(1,2) zeros(1,2);
%     -ones(1,2) zeros(1,2);
%     ones(1,2) zeros(1,2);
%     -ones(1,2) zeros(1,2)];

% Inicial conditions
f = zeros(2*N,1);
r = zeros(2*N,1);
y0 = zeros(2,1);
du0 = zeros(2,1);
u0 = zeros(2,1);
x0 = zeros(length(A), 1);
e0 = zeros(2,1);
[~,z01] = filter(b1(2:end),a1,0);
[~,z02] = filter(b2(2:end),a2,0);
[~,z03] = filter(b3(2:end),a3,0);
[~,z04] = filter(b4(2:end),a4,0);
% return
for k=1:qntd_k

    y(:,k) = y0;
    x(:,k) = A*x0 + B*du0 + D*e0;
    e(:,k) = y(:,k) - H*x(:,k);
    
    % free response
    for i=1:2*N
       f(i,1) = F(i,:)*x(:,k) + E(i,:)*e(:,k); 
    end

    % references
    if k>qntd_k/2
        r1 = 0.4;
    end
    r = repmat([r1; r2], ceil(2*N/2), 1);
    r = r(1:2*N);
    ref(:,k) = r;
    
    % quadratic equation bb
    bb = 2*(f-r)'*Q*G;
    
    % constraints matrix Bc
    % Bcon = [y_max-f;
    %     -y_min+f;
    %     (u_max-u0);
    %     -(u_min-u0);
    %     du_max;
    %     -du_min];

    % QP solver
    opt = optimoptions('quadprog','Display','off');
    Acon = []; Bcon = []; % -> let uncomment to remove constraints
    sol = quadprog(HH,bb,Acon,Bcon,[],[],[],[],[],opt);

    du(:,k) = sol(1:2);
    u(:,k) = u0 + du(:,k);

    [y1,z01] = filter(b1(2:end),a1,u(1,k),z01);
    [y2,z02] = filter(b2(2:end),a2,u(2,k),z02);
    Y1 = y1+y2;

    [y3,z03] = filter(b3(2:end),a3,u(1,k),z03);
    [y4,z04] = filter(b4(2:end),a4,u(2,k),z04);
    Y2 = y3+y4;

    du0 = du(:,k);
    u0 = u(:,k);
    x0 = x(:,k);
    e0 = e(:,k);
    y0 = [Y1; Y2];
end

% return
%% Figure
t = 1:1:qntd_k;

figure
subplot(2,1,1)
hold on
grid on
plot(t,ref(1,:), '--k', 'linewidth', 1)
plot(t,ref(2,:), '--k', 'linewidth', 1)
plot(t,y(1,:), 'r', 'linewidth', 1)
plot(t,y(2,:), 'b', 'linewidth', 1)
xlabel('Samples')
ylabel('y')
title('Output')

subplot(2,1,2)
hold on
grid on
plot(t,u(1,:), 'r', 'linewidth', 1)
plot(t,u(2,:), 'b', 'linewidth', 1)
xlabel('Samples')
ylabel('u')
title('Control signal')

sgtitle("\textbf{Output and Control Signals for N = }" + N + "\textbf{ and Nu = }" + Nu)

fprintf('\n\n Code Author: Jose Sergio Cruz Dantas Junior - UFC \n')
% fim