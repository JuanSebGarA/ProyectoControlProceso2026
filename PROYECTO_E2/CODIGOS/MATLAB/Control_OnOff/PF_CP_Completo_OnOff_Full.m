clc; clear; close all

%% Parametros modelo
% Parámetros maximos
pmax = [0.0418 0.2109 0.6995 0.1197 0.0762]; % mu,q,l,rho,pi
% Parámetros minimos
pmin = [0.0196 0.0006]; % q,l
% Parámetros modelo
param_mod = [0.9597 0.1908 0.0167 0.1002 0.579 12.5596 66.5337 100]; % Yxs, Yls, ms, mu_max, qmin, qmax, lmin, lmax, rho_max, pi_max, KS1, KN, KL, KI, KI2

% Flujos iniciales
F_0 = [0.25 0 0 0.025 0.025];

F_in = 2.9843; %L/h

% Entradas
x_in = 0;
S1_in = 60;
S2_in = 20;
Q_in = 0;
L_in = 0;
I0_in = 200;

% Otros parámetros no entradas
R = 0.1588;
alpha = 0.0254804;
B = 0.04;
V = 270;

param_in = [I0_in,R,x_in,S1_in,S2_in,Q_in,L_in,alpha,B,V]; %I0, R, x_in, S1_in, S2_in, Q_in, L_in, alpha, V

u = [F_in,x_in,S1_in,S2_in,Q_in,L_in,I0_in]'; %Igual que la de arriba pero para la S-funchon


tspan = linspace(0,200000,120); % Tiempo infinito

[t,y_out] = ode15s(@(t,Y) odeset(t,Y,pmax,pmin,param_in,param_mod,F_in),tspan,F_0);

x_out = y_out(:,1)';
S1_out = y_out(:,2)';
S2_out = y_out(:,3)';
Q_out = y_out(:,4)';
L_out = y_out(:,5)';
X_out = x_out + Q_out + L_out;

ODE_ss = [y_out(end,1),y_out(end,2),y_out(end,3),y_out(end,4),y_out(end,5)];

SS_desv = (F_0-ODE_ss)';

L_desv = (L_out-L_out(1))';

L_ss = ODE_ss(5);


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

t_step = 2000;

t_sim = 12000;

t_sp_change = 5000;

%% Flujo

%% SPT DR

% L_ss_ini = L_ss * 0.6;
% L_ss_fin = L_ss;
% 
% H_L = 2;
% 
% out = sim('PF_implementacion_OnOff_Flujo');
% 
% F_in_t = out.Inputs.Data(:,1);
% S1_in_t = out.Inputs.Data(:,3);
% 
% L_ODE = out.Outputs_Reactor.Data(:,5);
% 
% L_error = out.L_set.Data(:,1);
% 
% time_inputs = out.Inputs.Time;
% 
% L_ss_track = out.L_ss_track.Data(:,1);
% 
% 
% figure(1)
% 
% subplot(2,1,1)
% plot(time_inputs,L_ODE,'LineWidth',2.5)
% hold on
% 
% plot(time_inputs,L_ss_track,'r--','LineWidth',2.5)
% 
% plot(time_inputs,L_ss_track + H_L/2,'k--','LineWidth',1)
% plot(time_inputs,L_ss_track - H_L/2,'k--','LineWidth',1)
% 
% ylabel('L')
% title('Respuesta del sistema y banda de histéresis')
% legend('L(t)','L_{set}','L_{sup}','L_{inf}', Location='northwest')
% grid on
% 
% subplot(2,1,2)
% 
% yyaxis left
% stairs(time_inputs,F_in_t,'LineWidth',2.5)
% ylabel('F_{in} (L/h)')
% 
% yyaxis right
% stairs(time_inputs,S1_in_t,'LineWidth',2.5)
% ylabel('S_{1,in} (g/L)')
% 
% xlabel('Tiempo')
% title('Acción del controlador y perturbación')
% legend('F_{in}','S_{1,in}','Location','best')
% grid on
% set(gcf,'Color','w')

% exportgraphics(gcf,...
%     'Disturbance_Rejection_OnOff_Flujo_08.png',...
%     'Resolution',600)

% Video

% v = VideoWriter('Disturbance_Rejection_OnOff_Flujo.mp4','MPEG-4');
% v.FrameRate = 20;
% open(v)
% 
% fig = figure('Position',[100 100 1000 700]);
% 
% for k = 1:length(time_inputs)
% 
%     clf
% 
%     subplot(2,1,1)
% 
%     plot(time_inputs(1:k),L_ODE(1:k),'LineWidth',2.5)
%     hold on
% 
%     plot(time_inputs(1:k),L_ss_track(1:k),'r--','LineWidth',2.5)
% 
%     plot(time_inputs(1:k), ...
%         L_ss_track(1:k)+H_L/2,...
%         'k--','LineWidth',1)
% 
%     plot(time_inputs(1:k), ...
%         L_ss_track(1:k)-H_L/2,...
%         'k--','LineWidth',1)
% 
%     ylabel('L')
%     title('Respuesta del sistema y banda de histéresis')
%     legend('L(t)','L_{set}','L_{sup}','L_{inf}',...
%         'Location','northwest')
% 
%     grid on
% 
%     xlim([min(time_inputs) max(time_inputs)])
% 
%     subplot(2,1,2)
% 
%     yyaxis left
%     stairs(time_inputs(1:k),F_in_t(1:k),...
%         'LineWidth',2.5)
%     ylabel('F_{in} (L/h)')
% 
%     yyaxis right
%     stairs(time_inputs(1:k),S1_in_t(1:k),...
%         'LineWidth',2.5)
%     ylabel('S_{1,in} (g/L)')
% 
%     xlabel('Tiempo')
%     title('Acción del controlador y perturbación')
% 
%     grid on
% 
%     xlim([min(time_inputs) max(time_inputs)])
% 
%     frame = getframe(fig);
%     writeVideo(v,frame);
% 
% end
% 
% close(v)


%% Ruido Sensor NO SIRVE

% H_L_R = 0.2;
% 
% 
% out_R = sim('PF_implementacion_OnOff_Flujo_Ruido');
% 
% F_in_t_R = out_R.Inputs.Data(:,1);
% 
% L_ODE_R = out_R.Outputs_Reactor.Data(:,5);
% 
% L_error_R = out_R.L_set.Data(:,1);
% 
% time_inputs_R = out_R.Inputs.Time;
% 
% 
% 
% figure(1)
% 
% subplot(2,1,1)
% plot(time_inputs_R,L_ODE_R,'LineWidth',1.5)
% hold on
% 
% yline(L_ss,'r--','L_{set}','LineWidth',1.5)
% 
% yline(L_ss + H_L_R/2,'k--','L_{sup}','LineWidth',1)
% yline(L_ss - H_L_R/2,'k--','L_{inf}','LineWidth',1)
% 
% ylabel('L')
% title('Respuesta del sistema y banda de histéresis')
% legend('L(t)','L_{set}','L_{sup}','L_{inf}', Location='northwest')
% grid on
% 
% subplot(2,1,2)
% stairs(time_inputs_R,F_in_t_R,'LineWidth',1.5)
% ylabel('F_in')
% xlabel('Tiempo')
% title('Acción del controlador On-Off')
% grid on

%% S1

% L_ss_ini = L_ss * 0.6;
% L_ss_fin = L_ss;
% 
% H_L = 1;
% 
% out_S1 = sim('PF_implementacion_OnOff_S1');
% 
% S1_in_t_S1 = out_S1.Inputs.Data(:,3);
% I0_in_t_S1 = out_S1.Inputs.Data(:,7);
% 
% L_ODE_S1 = out_S1.Outputs_Reactor.Data(:,5);
% 
% L_error_S1 = out_S1.L_set.Data(:,1);
% 
% time_inputs_S1 = out_S1.Inputs.Time;
% 
% L_ss_track_S1 = out_S1.L_ss_track.Data(:,1);
% % 
% % 
% figure(1)
% 
% subplot(2,1,1)
% plot(time_inputs_S1,L_ODE_S1,'LineWidth',2.5)
% hold on
% 
% plot(time_inputs_S1,L_ss_track_S1,'r--','LineWidth',2.5)
% 
% plot(time_inputs_S1,L_ss_track_S1 + H_L/2,'k--','LineWidth',1)
% plot(time_inputs_S1,L_ss_track_S1 - H_L/2,'k--','LineWidth',1)
% 
% ylabel('L')
% title('Respuesta del sistema y banda de histéresis')
% legend('L(t)','L_{set}','L_{sup}','L_{inf}', Location='southeast')
% grid on
% 
% subplot(2,1,2)
% 
% yyaxis left
% stairs(time_inputs_S1,S1_in_t_S1,'LineWidth',2.5)
% ylabel('S_{1,in} (g/L)')
% 
% yyaxis right
% stairs(time_inputs_S1,I0_in_t_S1,'LineWidth',2.5)
% ylabel('I_{0}')
% 
% xlabel('Tiempo')
% title('Acción del controlador y perturbación')
% legend('S_{1,in}','I_0','Location','best')
% grid on
% set(gcf,'Color','w')
% 
% exportgraphics(gcf,...
%     'Disturbance_Rejection_OnOff_S1_Corregida.png',...
%     'Resolution',600)

% Video

% v = VideoWriter('Disturbance_Rejection_OnOff_S1.mp4','MPEG-4');
% v.FrameRate = 40;
% open(v)
% 
% fig = figure('Position',[100 100 1000 700]);
% 
% for k = 1:length(time_inputs_S1)
% 
%     clf
% 
%     subplot(2,1,1)
% 
%     plot(time_inputs_S1(1:k),L_ODE_S1(1:k),'LineWidth',2.5)
%     hold on
% 
%     plot(time_inputs_S1(1:k),L_ss_track_S1(1:k),'r--','LineWidth',2.5)
% 
%     plot(time_inputs_S1(1:k), ...
%         L_ss_track_S1(1:k)+H_L/2,...
%         'k--','LineWidth',1)
% 
%     plot(time_inputs_S1(1:k), ...
%         L_ss_track_S1(1:k)-H_L/2,...
%         'k--','LineWidth',1)
% 
%     ylabel('L')
%     title('Respuesta del sistema y banda de histéresis')
%     legend('L(t)','L_{set}','L_{sup}','L_{inf}',...
%         'Location','northwest')
% 
%     grid on
% 
%     xlim([min(time_inputs_S1) max(time_inputs_S1)])
% 
%     subplot(2,1,2)
% 
%     yyaxis left
%     stairs(time_inputs_S1(1:k),S1_in_t_S1(1:k),...
%         'LineWidth',2.5)
%     ylabel('S1_{in} (g/L)')
% 
%     yyaxis right
%     stairs(time_inputs_S1(1:k),I0_in_t_S1(1:k),...
%         'LineWidth',2.5)
%     ylabel('I_{0}')
% 
%     xlabel('Tiempo')
%     title('Acción del controlador y perturbación')
% 
%     grid on
% 
%     xlim([min(time_inputs_S1) max(time_inputs_S1)])
% 
%     frame = getframe(fig);
%     writeVideo(v,frame);
% 
% end
% 
% close(v)



%% Intensidad

%% SPT and DR

%SPT

L_ss_ini = L_ss * 0.6;
L_ss_fin = L_ss;

H_L_I = 2;

out_I = sim('PF_implementacion_OnOff_Intensidad');

I0_in_t_I = out_I.Inputs.Data(:,7);
S1_in_t_I = out_I.Inputs.Data(:,3);

L_ODE_I = out_I.Outputs_Reactor.Data(:,5);

L_error_I = out_I.L_set.Data(:,1);

time_inputs_I = out_I.Inputs.Time;

L_ss_track_I = out_I.L_ss_track_I.Data(:,1);

figure(1)

subplot(2,1,1)
plot(time_inputs_I,L_ODE_I,'LineWidth',2.5)
hold on

plot(time_inputs_I,L_ss_track_I,'r--','LineWidth',2.5)

plot(time_inputs_I,L_ss_track_I + H_L_I/2,'k--','LineWidth',1)
plot(time_inputs_I,L_ss_track_I - H_L_I/2,'k--','LineWidth',1)

ylabel('L')
title('Respuesta del sistema y banda de histéresis')
legend('L(t)','L_{set}','L_{sup}','L_{inf}', Location='northwest')
grid on

subplot(2,1,2)

yyaxis left
stairs(time_inputs_I,I0_in_t_I,'LineWidth',2.5)
ylabel('I0_{in}')

yyaxis right
stairs(time_inputs_I,S1_in_t_I,'LineWidth',2.5)
ylabel('S_{1,in} (g/L)')

xlabel('Tiempo')
title('Acción del controlador y perturbación')
legend('I0_{in}','S_{1,in}','Location','best')
grid on
set(gcf,'Color','w')


% exportgraphics(gcf,...
%     'Disturbance_Rejection_OnOff_Intensidad_Corregido.png',...
%     'Resolution',600)

% Video
% v = VideoWriter('Disturbance_Rejection_OnOff_Intensidad.mp4','MPEG-4');
% v.FrameRate = 20;
% open(v)
% 
% fig = figure('Position',[100 100 1000 700]);
% 
% for k = 1:length(time_inputs_I)
% 
%     clf
% 
%     subplot(2,1,1)
% 
%     plot(time_inputs_I(1:k),L_ODE_I(1:k),'LineWidth',2.5)
%     hold on
% 
%     plot(time_inputs_I(1:k),L_ss_track_I(1:k),'r--','LineWidth',2.5)
% 
%     plot(time_inputs_I(1:k), ...
%         L_ss_track_I(1:k)+H_L_I/2,...
%         'k--','LineWidth',1)
% 
%     plot(time_inputs_I(1:k), ...
%         L_ss_track_I(1:k)-H_L_I/2,...
%         'k--','LineWidth',1)
% 
%     ylabel('L')
%     title('Respuesta del sistema y banda de histéresis')
%     legend('L(t)','L_{set}','L_{sup}','L_{inf}',...
%         'Location','northwest')
% 
%     grid on
% 
%     xlim([min(time_inputs_I) max(time_inputs_I)])
% 
%     subplot(2,1,2)
% 
%     yyaxis left
%     stairs(time_inputs_I(1:k),I0_in_t_I(1:k),...
%         'LineWidth',2.5)
%     ylabel('I_{0}')
% 
%     yyaxis right
%     stairs(time_inputs_I(1:k),S1_in_t_I(1:k),...
%         'LineWidth',2.5)
%     ylabel('S_{1,in} (g/L)')
% 
%     xlabel('Tiempo')
%     title('Acción del controlador y perturbación')
% 
%     grid on
% 
%     xlim([min(time_inputs_I) max(time_inputs_I)])
% 
%     frame = getframe(fig);
%     writeVideo(v,frame);
% 
% end
% 
% close(v)

%% Ruido Sensor

% H_L_IR = 1;
% 
% 
% out_IR = sim('PF_implementacion_OnOff_Intensidad_Ruido');
% 
% I0_in_t_IR = out_IR.Inputs.Data(:,7);
% 
% L_ODE_IR = out_IR.Outputs_Reactor.Data(:,5);
% 
% L_error_IR = out_IR.L_set.Data(:,1);
% 
% time_inputs_IR = out_IR.Inputs.Time;
% 
% 
% 
% figure(1)
% 
% subplot(2,1,1)
% plot(time_inputs_IR,L_ODE_IR,'LineWidth',1.5)
% hold on
% 
% yline(L_ss,'r--','L_{set}','LineWidth',1.5)
% 
% yline(L_ss + H_L_IR/2,'k--','L_{sup}','LineWidth',1)
% yline(L_ss - H_L_IR/2,'k--','L_{inf}','LineWidth',1)
% 
% ylabel('L')
% title('Respuesta del sistema y banda de histéresis')
% legend('L(t)','L_{set}','L_{sup}','L_{inf}', Location='northwest')
% grid on
% 
% subplot(2,1,2)
% stairs(time_inputs_IR,I0_in_t_IR,'LineWidth',1.5)
% ylabel('F_in')
% xlabel('Tiempo')
% title('Acción del controlador On-Off')
% grid on







%%

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

