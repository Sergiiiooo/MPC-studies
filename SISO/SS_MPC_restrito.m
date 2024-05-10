clc
clear
close all
warning off

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex');

Ts = 1;
s = tf('s');
z = tf('z', Ts);

Gz = 0.2/(1 - 0.9*z^-1);
% Gz=c2d(1/(1*s+1)^2,Ts,'zoh');


Delta = z/(z-1);

[numGz, denGz] = tfdata(Gz, 'v');
[numGzD, denGzD] = tfdata(Gz*Delta, 'v');

% c polynomial
alpha = 0;
beta = 0;
c = [1 -alpha beta];
c = [c zeros(1,length(denGzD)-length(c))];

%% State space
[A,H,B,~] = tf2ss(numGz,denGzD);
B = B';
A = A';
H = H';

[~,~,D,~] = tf2ss(c,denGzD);
D = D';

%% Controle

% Horizonte de predicao
N = 3;

% Horizonte de controle
Nu = 3;

% G
for i=1:N
   G(i,1) = H*A^(i-1)*B;
end

for i=2:N
   for j=2:Nu
      G(i,j) = G(i-1,j-1);
   end
end

% F
for i=1:N
   F(i,:) = H*A^(i); 
end

% E 
for i=1:N
   E(i,1) = H*A^(i-1)*D; 
end

%% Simulacao
% parametros
Tsim = 100;
Tref = 0;
Tdist = 50;
Tnoise = 0.9*Tsim;
RefAmplitude = 1;
RefInitialValue = 0;
DistAmplitude = -0.25;
NoisePower = 1e-2;

ref = RefAmplitude*ones(ceil(Tsim/Ts),1);
q = [zeros(ceil(Tdist/Ts),1); DistAmplitude*ones(ceil((Tsim-Tdist)/Ts),1)];

fu = zeros(N,1);
f = zeros(N,1);
r = zeros(N,1);

yf = 0;
y0 = 0;
du0 = 0;
u0 = 0;
x0 = zeros(length(A), 1);
e0 = 0;
[~,z0] = filter(numGz,denGz,0);
b0 = [1;zeros(Nu-1,1)];

for k=1:ceil(Tsim/Ts)
    y(k) = y0;
    
    x(:,k) = A*x0 + B*du0 + D*e0; %x(k) = Ax(k-1) + Bdu(k-1) + De(k-1)
    e(k) = y(k) - H*x(:,k); %e(k)=yk)-\hat{y}(k)
    
    % resposata livre
    for i=1:N
       f(i,1) = F(i,:)*x(:,k) + E(i,1)*e(k); 
    end
    
    % ref
    for i=1:N
        r(i,1) = RefAmplitude; % case the reference is known
    end
    
    % pesos de Qdelta e Qlambda
    delta_i = 1;
    lambda_i = 1;

    % Qdelta
    Qdelta = delta_i*eye(N);

    % Qlambda
    Qlambda = lambda_i*eye(Nu);
    
    % equacao quadrtica HH
    HH = 2*G'*Qdelta*G + 2*Qlambda;
    
    % equacao quadrtica b
    b = 2*G'*Qdelta*(f-r);

    % constraints
    ymax = 1.2;
    ymin = 0;
    u_max = 0.9*ones(Nu,1);
    u_min = 0*ones(Nu,1);
    du_max = 1*ones(Nu,1);
    du_min = -1*ones(Nu,1);
    T = tril(ones(Nu));
    
    % constraints' matrices
    Acon1 = [G;-G];
    Acon2 = [T;-T];
    Acon3 = [eye(Nu); -eye(Nu)];
    Bcon1 = [ymax-f; f-ymin];
    Bcon2 = [u_max-T*b0*u0; T*b0*u0-u_min];
    Bcon3 = [du_max; -du_min];

    Acon = [Acon1; Acon2; Acon3];
    Bcon = [Bcon1; Bcon2; Bcon3];
    
    % resolvedor
    opt = optimoptions('quadprog','Display', 'off');
    sol = quadprog(HH,b,Acon,Bcon,[],[],[],[],[],opt);
    % sol = quadprog(HH,b,[],[],[],[],[],[],[],opt);

K = (G'*Qdelta*G + Qlambda)\G'*Qdelta;
K = K(1,:);

% Reference gain Kr
Kr = sum(K);


   
    du(k) = sol(1,1);
    % du(k) = Kr*r-K*F(i,:)*x(:,k) - K*E(i,1)*e(k);
    u(k) = u0 + du(k);

    ut = u(k) + q(k);
    [yf,zf] = filter(numGz,denGz,ut,z0);

    y0 = yf;
    z0 = zf;
    du0 = du(k);
    u0 = u(k);
    x0 = x(:,k);
    e0 = e(k);
    
end

figure
subplot(3,1,1)
plot(ref, '--k', 'linewidth', 2)
hold on
stairs(y, 'b', 'linewidth', 1)
grid on
axis tight
legend('Referencia, r(k)', 'saida, y(t)', 'location', 'best')
subplot(3,1,2)
stairs(u, 'b', 'linewidth', 1)
grid on
axis tight
legend('Sinal de controle, u(k)', 'location', 'best')
subplot(3,1,3)
stairs(du, 'b', 'linewidth', 1)
grid on
axis tight
xlabel('Samples, k')
legend('Variacao de controle, $\Delta$u(k)', 'location', 'best')

figure
subplot(2,1,1)
stairs(x(:,1), 'b', 'linewidth', 1)
hold on
stairs(x(:,2), 'r', 'linewidth', 1)
grid on
axis tight
legend('$x_{1}$(k)', '$x_{2}$(k)',  'location', 'best')
subplot(2,1,2)
stairs(e, 'b', 'linewidth', 1)
grid on
axis tight
xlabel('Samples, k')
legend('e(k)', 'location', 'best')

return


figure
subplot(3,1,1)
plot(ref, '--k', 'linewidth', 2)
hold on
stairs(y, 'b', 'linewidth', 1)
stairs(out.y, 'magenta', 'linewidth', 1)
grid on
axis tight
legend('Referencia, r(k)', 'saida, y(t)', 'location', 'best')
subplot(3,1,2)
stairs(u, 'b', 'linewidth', 1)
stairs(out.u, 'magenta', 'linewidth', 1)
grid on
axis tight
legend('Sinal de controle, u(k)', 'location', 'best')
subplot(3,1,3)
stairs(du, 'b', 'linewidth', 1)
stairs(out.du, 'magenta', 'linewidth', 1)
grid on
axis tight
xlabel('Samples, k')
legend('Variacao de controle, $\Delta$u(k)', 'location', 'best')

