clc; clear; close all

%% ================== PARÁMETROS ==================
param_max = [0.0418 0.2109 0.6995 0.1197 0.0762]; % mu,q,l,rho,pi
param_min = [0.0196 0.0006]; % q,l

param_in = [200,0.15,0,60,20,0,0,0.0254804,0.04,270]; % %I0, R, x_in, S1_in, S2_in, Q_in, L_in, alpha, V
param_mod = [0.9597 0.1908 0.0167 0.1002 0.579 12.5596 66.5337 100];

F_0 = [0.25 0 0 0.025 0.025];

tspan = linspace(0,1250,150);

F_in_opt = 2; % puedes cambiarlo o hacerlo vector después

%% ================== OPTIMIZACIÓN ==================
lb = [50, 0.5];   % límites
ub = [500, 3];

p0 = [250,0.15];

[p_opt, fval] = fmincon(@(p) objetivo_L(p,param_max,param_min,param_in,param_mod,F_0,tspan,F_in_opt),...
                        p0,[],[],[],[],lb,ub);

I0_opt = p_opt(1)
R_opt  = p_opt(2)
L_max  = -fval

%% ================== SIMULACIÓN ÓPTIMA ==================
param_in(1) = I0_opt;
param_in(2) = R_opt;

[t,y] = ode15s(@(t,Y) modelo_ode(t,Y,param_max,param_min,param_in,param_mod,F_in_opt),tspan,F_0);

x = y(:,1);
Q = y(:,4);
L = y(:,5);
X = x + Q + L;

%% ================== GRÁFICA ==================
figure
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1])

plot(t,X,'k','LineWidth',2); hold on
plot(t,L,'g','LineWidth',2)
plot(t,Q,'r--','LineWidth',2)
plot(t,x,'b','LineWidth',2)

legend('X_{total}','L','Q','x','Location','best')
xlabel('Tiempo [h]','FontSize',12)
ylabel('Concentración [g/L]','FontSize',12)
title(sprintf('Óptimo: I0 = %.2f | R = %.2f',I0_opt,R_opt))
set(gca,'FontSize',12)

%% ================== FUNCIÓN OBJETIVO ==================
function J = objetivo_L(p, param_max, param_min, param_in, param_mod, F_0, tspan, F_in)

I0 = p(1);
R  = p(2);

param_in(1) = I0;
param_in(2) = R;

[t,y] = ode15s(@(t,Y) modelo_ode(t,Y,param_max,param_min,param_in,param_mod,F_in),tspan,F_0);

L = y(:,5);

J = -max(L); % maximizar lípidos

end

%% ================== MODELO ==================
function F = odeset(t,Y,max,min,param_in,param_mod,F_in)

% Entradas
x = Y(1); S1 = Y(2); S2 = Y(3); Q = Y(4); L = Y(5);

%Parametros Max-Min
mu_max = max(1); q_max = max(2); l_max = max(3); rho_max = max(4); pi_max = max(5);
q_min = min(1); l_min = min(2);

% Parametros Entrada
I0 = param_in(1); R = param_in(2); x_in = param_in(3); S1_in = param_in(4); S2_in = param_in(5);
Q_in = param_in(6); L_in = param_in(7); alpha = param_in(8); B = param_in(9); V = param_in(10);

% Parametros modelos
Yxs = param_mod(1); Yls = param_mod(2); ms = param_mod(3); K_S1 = param_mod(4);
K_N = param_mod(5); K_L = param_mod(6); K_I = param_mod(7); K_I2 = param_mod(8);

% Totales
X = x + Q + L;
q = Q/X;
l = L/X;

%I = (((2*I0)/(alpha*B*X*1000))*(1/R^2))*(R-((1)/(alpha*B*X*1000))+((exp(-alpha*B*X*R*1000))/(alpha*B*X*1000)));
I = (((2*I0)/(alpha*B*X))*(1/R^2))*(R-((1)/(alpha*B*X))+((exp(-alpha*B*X*R))/(alpha*B*X)));

% Parametros ODE
mu = mu_max * (1 - (q_min/q)) * (1 - (l_min/l)) * (S1/(K_S1+S1)) * (I/(K_I + I + I^2/K_I2));
mu = mu * (mu >= 0) * (q >= 0) * (l >= 0) * (q >= q_min) * (l >= l_min) * (S1 >= 0);

rho = rho_max * (S2/(K_N + S2)) * ((q_max-q)/(q_max - q_min));
rho = rho *( rho >=0) *( S2 >=0) *( q >=0) *( S2 >=0);

pi = pi_max * (S1/(K_L+S1)) * (1 - q) * ((l_max - l)/l_max);
pi = pi *( pi >=0) *( S1 >=0) *( q >=0) *( l >=0) *( q < q_max ) *( l < l_max );

% ODE
dx_dt = mu * x + (F_in/V)*(x_in - x);
dS1_dt = -(1/Yxs)*(mu * x) - (1/Yls)*(pi*x) - ms*x + (F_in/V)*(S1_in-S1);
dS2_dt = - rho*x + (F_in/V)*(S2_in-S2);
dQ_dt = rho*x - mu * Q + (F_in/V) * (Q_in - Q);
dL_dt = pi * x - mu * L + (F_in/V) * (L_in - L);


F = [dx_dt, dS1_dt, dS2_dt, dQ_dt, dL_dt]';
end
