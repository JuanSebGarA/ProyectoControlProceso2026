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

[A1, B1, C1, D1] = linmod('PF_implementacion_modeloL2',ODE_ss,u);

SS_desv = (F_0-ODE_ss)';

L_desv = (L_out-L_out(1))';

%% Análisis de respuestas
Sys_ss = ss(A1,B1,C1,D1);   %Objeto-espacio de estados % D,'StateName',{'Ca';'Cb'},'InputName',{'D';'Caf'},'OutputName',{'Ca';'Cb'}
Sys_tf = tf(Sys_ss);        %Transformación a función de tranferencia

%% Polos

den = Sys_tf.den{1};
Polos = roots(den);
Polos2 = pole(Sys_ss);


error_polos = abs(Polos - Polos2)/100;

%% Estabilidad
v_prop = eig(A1);
cont = ctrb(A1,B1);

controlabilidad = length(A1) - rank(cont);

obs = obsv(A1,C1);
observabilidad = length(A1) - rank(obs);

%% Reducción

Sys_min = minreal(Sys_ss);

mSys_min = sminreal(Sys_ss);

rSys_min = balred(Sys_tf,2);

[brSys_min, g]= balreal(Sys_ss);

elim = (g<0.3); % Identificación de VSH menores a 0.3
modSys_min = modred(brSys_min,elim,'MatchDC');
%% RGA

for i = 1:length(B1)
    for j = 1:length(A1)
        K1p(i,j) = - C1(j,:)*inv(A1)*B1(:,i)+D1(j,i);
    end
end


for j=1:length(B1)
    [num1,den1]=ss2tf(A1,B1,C1,D1,j);
    for i=1:length(A1)
        K_f(j,i)=poly2sym(num1(i,:))/poly2sym(den1);
    end
end
x = 0;
K1 = eval(K_f);

RGA = K1.*(pinv(K1))';
RGA1 = K1p.*(pinv(K1p))';

%% SVD

[U, S, V1] = svd(K1);
CN = S(1,1)/S(5,5);
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

% Entradas
F_in_t = out.Inputs.Data(:,1);
x_in_t = out.Inputs.Data(:,2);
S1_in_t = out.Inputs.Data(:,3);
S2_in_t = out.Inputs.Data(:,4);
Q_in_t = out.Inputs.Data(:,5);
L_in_t = out.Inputs.Data(:,6);
I0_in_t = out.Inputs.Data(:,7);

time_inputs = out.Inputs.Time;

% ODES

x_ODE = out.Outputs_Reactor.Data(:,1);
S1_ODE = out.Outputs_Reactor.Data(:,2);
S2_ODE = out.Outputs_Reactor.Data(:,3);
Q_ODE = out.Outputs_Reactor.Data(:,4);
L_ODE = out.Outputs_Reactor.Data(:,5);
time_ODES = out.Outputs_Reactor.Time;

% Modelo lineal

x_Lin = out.lineal_model.Data(:,1);
S1_Lin = out.lineal_model.Data(:,2);
S2_Lin = out.lineal_model.Data(:,3);
Q_Lin = out.lineal_model.Data(:,4);
L_Lin = out.lineal_model.Data(:,5);
time_Lin = out.lineal_model.Time;


Dist_S1 = S1_in_t-S1_in;


%% Figuras

%%% FIGURA 1

figure(1)
set(gcf,'Units','centimeters','Position',[0 0 16 25])

% Configuración global
set(groot, 'defaultAxesFontName', 'Palatino Linotype')
set(groot, 'defaultAxesFontSize', 11)
set(groot, 'defaultLineLineWidth', 2.5)

%----------------------
subplot(5,1,1)
h1 = plot(time_ODES, x_ODE,'k'); hold on
h2 = plot(time_Lin, x_Lin, 'r--');
ylabel('x [g/L]')

%----------------------
subplot(5,1,2)
plot(time_ODES, S1_ODE,'k'); hold on
plot(time_Lin, S1_Lin, 'r--')
ylabel('S_1 [g/L]')

%----------------------
subplot(5,1,3)
plot(time_ODES, S2_ODE,'k'); hold on
plot(time_Lin, S2_Lin, 'r--')
ylabel('S_2 [g/L]')

%----------------------
subplot(5,1,4)
plot(time_ODES, Q_ODE,'k'); hold on
plot(time_Lin, Q_Lin, 'r--')
ylabel('Q [g/L]')

%----------------------
subplot(5,1,5)
plot(time_ODES, L_ODE,'k'); hold on
plot(time_Lin, L_Lin, 'r--')
ylabel('L [g/L]')
xlabel('Tiempo [h]')


lgd = legend([h1 h2], {'Modelo no lineal','Modelo linealizado'}, ...
    'Orientation','horizontal');

lgd.FontName = 'Palatino Linotype';
lgd.Units = 'normalized';
lgd.Position = [0.3 0.01 0.4 0.04]; % abajo centrada


% sgtitle('Perturbación: S_2 x1.5 ','FontName','Palatino Linotype','FontSize',12)
% %print(gcf,'perturbacion_3','-dpng','-r300')

%%% FIGURA 2

% figure(2)
% set(gcf,'Units','centimeters','Position',[2 2 12 8]) % tamaño compacto
%
% % Configuración
% set(groot, 'defaultAxesFontName', 'Palatino Linotype')
% set(groot, 'defaultAxesFontSize', 11)
% set(groot, 'defaultLineLineWidth', 2.5)
%
% % Graficas
% h1 = plot(time_ODES, L_ODE,'g'); hold on
% h2 = plot(time_Lin, L_Lin,'g--');
%
% h3 = plot(time_ODES, Q_ODE,'r');
% h4 = plot(time_Lin, Q_Lin,'r--');
%
% h5 = plot(time_ODES, x_ODE,'b');
% h6 = plot(time_Lin, x_Lin,'b--');
%
% % Labels
% xlabel('Tiempo [h]','FontName','Palatino Linotype')
% ylabel('Concentración [g/L]','FontName','Palatino Linotype')
%
% % Leyenda clara
% legend([h1 h2 h3 h4 h5 h6], ...
%     {'L (no lineal)','L (lineal)', ...
%      'Q (no lineal)','Q (lineal)', ...
%      'x (no lineal)','x (lineal)'}, ...
%      'Location','eastoutside', ...
%      'FontName','Palatino Linotype','FontSize',10)
%
% set(gca,'FontName','Palatino Linotype')
%
% % Exportar
% print(gcf,'comparacion_variables','-dpng','-r300')

% figure(3)
% set(groot, 'defaultLineLineWidth', 0.5)
% pzmap(Sys_tf)
% grid on
% ylabel('Eje Imaginario','FontWeight','bold','FontSize',14)
% xlabel('Eje Real','FontWeight','bold','FontSize',14)
% set(gca,'FontSize',12,'FontWeight','bold')
% print(gcf,'pzmap','-dpng','-r300')
%


% figure(4)
% iopzmap(Sys_tf)
% ylabel('Eje Imaginario','FontWeight','bold','FontSize',14)
% xlabel('Eje Real','FontWeight','bold','FontSize',14)
% set(gca,'FontSize',12,'FontWeight','bold')

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

