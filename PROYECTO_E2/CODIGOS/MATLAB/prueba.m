clc; clear; close all

%% Parámetros
param_max = [0.0418 0.2109 0.6995 0.1197 0.0762];
param_min = [0.0196 0.0006];

param_in = [200,1,0,60,20,0,0,0.0254804,0.04,270];
param_mod = [0.9597 0.1908 0.0167 0.1002 0.579 12.5596 66.5337 100];

F_0 = [0.25 0.1 0.1 0.025 0.025];

%% Optimización
lb = [50 0.2 0];
ub = [500 5 10];

p0 = [200 1 2];

options_fmin = optimoptions('fmincon',...
                            'Display','iter',...
                            'Algorithm','sqp',...
                            'MaxFunctionEvaluations',5000,...
                            'MaxIterations',200);

[p_opt,fval] = fmincon(@(p) objetivo_L_ss(p,param_max,param_min,param_in,param_mod,F_0),...
                       p0,[],[],[],[],lb,ub,[],options_fmin);

I0_opt = p_opt(1)
R_opt  = p_opt(2)
F_opt  = p_opt(3)

L_max = -fval

%% Estado estable óptimo
param_in(1) = I0_opt;
param_in(2) = R_opt;

[~,y] = ode15s(@(t,Y) modelo_ode(t,Y,param_max,param_min,param_in,param_mod,F_opt),...
               [0 5000],F_0);

Yguess = y(end,:);

options_fsolve = optimoptions('fsolve',...
                              'Display','iter',...
                              'FunctionTolerance',1e-8,...
                              'StepTolerance',1e-8);

Yss = fsolve(@(Y) modelo_ss(Y,param_max,param_min,param_in,param_mod,F_opt),...
             Yguess,options_fsolve);

x_ss  = Yss(1)
S1_ss = Yss(2)
S2_ss = Yss(3)
Q_ss  = Yss(4)
L_ss  = Yss(5)

%% Simulación dinámica final
tspan = linspace(0,5000,300);

[t,y] = ode15s(@(t,Y) modelo_ode(t,Y,param_max,param_min,param_in,param_mod,F_opt),...
               tspan,F_0);

x = y(:,1);
Q = y(:,4);
L = y(:,5);

X = x + Q + L;

%% Gráfica
figure
set(gcf,'Units','centimeters','Position',[2 2 12 8])

plot(t,X,'k','LineWidth',2.5); hold on
plot(t,L,'g','LineWidth',2.5)
plot(t,Q,'r--','LineWidth',2.5)
plot(t,x,'b','LineWidth',2.5)

yline(L_ss,'r:','LineWidth',2.5)

grid on

legend('X_{total}','L','Q','x','Location','best','FontName','Palatino Linotype')

xlabel('Tiempo [h]','FontSize',11,'FontName','Palatino Linotype')
ylabel('Concentración [g/L]','FontSize',11,'FontName','Palatino Linotype')

set(gca,'FontSize',11,'FontName','Palatino Linotype')

%% Función objetivo
function J = objetivo_L_ss(p,param_max,param_min,param_in,param_mod,F_0)

I0 = p(1);
R = p(2);
F_in = p(3);

param_in(1) = I0;
param_in(2) = R;

try

    warning off

    [~,y] = ode15s(@(t,Y) modelo_ode(t,Y,param_max,param_min,param_in,param_mod,F_in),...
                   [0 5000],F_0);

    Yguess = y(end,:);

    options = optimoptions('fsolve',...
                           'Display','off',...
                           'FunctionTolerance',1e-8,...
                           'StepTolerance',1e-8);

    Yss = fsolve(@(Y) modelo_ss(Y,param_max,param_min,param_in,param_mod,F_in),...
                 Yguess,options);

    Lss = Yss(5);

    if isnan(Lss) || isinf(Lss) || Lss < 0

        J = 1e6;

    else

        J = -Lss;

    end

catch ME

    disp('ERROR EN OPTIMIZACION:')
    disp(ME.message)

    J = 1e6;

end

end

%% Modelo estado estacionario
function F = modelo_ss(Y,param_max,param_min,param_in,param_mod,F_in)

F = modelo_ode(0,Y,param_max,param_min,param_in,param_mod,F_in);

end

%% Modelo dinámico
function F = modelo_ode(~,Y,param_max,param_min,param_in,param_mod,F_in)

%% Estados protegidos
x  = max(Y(1),1e-8);
S1 = max(Y(2),1e-8);
S2 = max(Y(3),1e-8);
Q  = max(Y(4),1e-8);
L  = max(Y(5),1e-8);

%% Parámetros max-min
mu_max  = param_max(1);
q_max   = param_max(2);
l_max   = param_max(3);
rho_max = param_max(4);
pi_max  = param_max(5);

q_min = param_min(1);
l_min = param_min(2);

%% Parámetros entrada
I0    = param_in(1);
R     = param_in(2);

x_in  = param_in(3);
S1_in = param_in(4);
S2_in = param_in(5);

Q_in  = param_in(6);
L_in  = param_in(7);

alpha = param_in(8);
B     = param_in(9);
V     = param_in(10);

%% Parámetros modelo
Yxs = param_mod(1);
Yls = param_mod(2);
ms  = param_mod(3);

K_S1 = param_mod(4);
K_N  = param_mod(5);
K_L  = param_mod(6);

K_I  = param_mod(7);
K_I2 = param_mod(8);

%% Totales protegidos
X = max(x + Q + L,1e-8);

q = max(Q/X,1e-8);
l = max(L/X,1e-8);

%% Intensidad lumínica
I = (((2*I0)/(alpha*B*X))*(1/R^2))*...
    (R-(1/(alpha*B*X))+((exp(-alpha*B*X*R))/(alpha*B*X)));

I = max(I,1e-8);

%% Cinéticas
mu = mu_max*(1-(q_min/q))*(1-(l_min/l))*...
     (S1/(K_S1+S1))*...
     (I/(K_I+I+(I^2/K_I2)));

rho = rho_max*(S2/(K_N+S2))*...
       ((q_max-q)/(q_max-q_min));

pi = pi_max*(S1/(K_L+S1))*...
     (1-q)*...
     ((l_max-l)/l_max);

%% Protección
mu  = max(mu,0);
rho = max(rho,0);
pi  = max(pi,0);

%% EDOs
dx_dt = mu*x + (F_in/V)*(x_in-x);

dS1_dt = -(1/Yxs)*(mu*x) ...
         -(1/Yls)*(pi*x) ...
         - ms*x ...
         + (F_in/V)*(S1_in-S1);

dS2_dt = -rho*x ...
         + (F_in/V)*(S2_in-S2);

dQ_dt = rho*x ...
        - mu*Q ...
        + (F_in/V)*(Q_in-Q);

dL_dt = pi*x ...
        - mu*L ...
        + (F_in/V)*(L_in-L);

F = [dx_dt dS1_dt dS2_dt dQ_dt dL_dt]';

end