clc; clear; close all

%% Parametros modelo
% Parámetros maximos
max = [0.0418 0.2109 0.6995 0.1197 0.0762]; % mu,q,l,rho,pi
% Parámetros minimos
min = [0.0196 0.0006]; % q,l
% Parámetros modelo
param_mod = [0.9597 0.1908 0.0167 0.1002 0.579 12.5596 66.5337 100]; % Yxs, Yls, ms, mu_max, qmin, qmax, lmin, lmax, rho_max, pi_max, KS1, KN, KL, KI, KI2

% Flujos iniciales
F_0 = [0.25 0 0 0.025 0.025];

F_in = 2.9922; %L/h

% Entradas
x_in = 0;
S1_in = 60;
S2_in = 20;
Q_in = 0;
L_in = 0;
I0_in = 200;

% Otros parámetros no entradas
R = 2.7312;
alpha = 0.0254804;
B = 0.04;
V = 270;

param_in = [I0_in,R,x_in,S1_in,S2_in,Q_in,L_in,alpha,B,V]; %I0, R, x_in, S1_in, S2_in, Q_in, L_in, alpha, V

u = [F_in,x_in,S1_in,S2_in,Q_in,L_in,I0_in]'; %Igual que la de arriba pero para la S-funchon


tspan = linspace(0,20000,120); % Tiempo infinito

[t,y_out] = ode15s(@(t,Y) odeset(t,Y,max,min,param_in,param_mod,F_in),tspan,F_0);

x_out = y_out(:,1)';
S1_out = y_out(:,2)';
S2_out = y_out(:,3)';
Q_out = y_out(:,4)';
L_out = y_out(:,5)';
X_out = x_out + Q_out + L_out;

ODE_ss = [y_out(end,1),y_out(end,2),y_out(end,3),y_out(end,4),y_out(end,5)];

SS_desv = (F_0-ODE_ss)';

L_desv = (L_out-L_out(1))';


%% Entradas simulink

F_in_ini = F_in;
F_in_fin = F_in;

x_in_ini = x_in;
x_in_fin = x_in;

S1_in_ini = S1_in;
S1_in_fin = S1_in*1.5;


S2_in_ini = S2_in;
S2_in_fin = S2_in;

Q_in_ini = Q_in;
Q_in_fin = Q_in;

L_in_ini = L_in;
L_in_fin = L_in;

I0_in_ini = I0_in;
I0_in_fin = I0_in;

t_step = 10000;

t_sim = 20000;

out = sim('PF_implementacion_modeloNL2');







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

