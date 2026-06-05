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


tspan = linspace(0,20000,120); % Tiempo infinito

[t,y_out] = ode15s(@(t,Y) odeset(t,Y,pmax,pmin,param_in,param_mod,F_in),tspan,F_0);

x_out = y_out(:,1)';
S1_out = y_out(:,2)';
S2_out = y_out(:,3)';
Q_out = y_out(:,4)';
L_out = y_out(:,5)';
X_out = x_out + Q_out + L_out;

ODE_ss = [y_out(end,1),y_out(end,2),y_out(end,3),y_out(end,4),y_out(end,5)];

L_ss = ODE_ss(5);

[A1, B1, C1, D1] = linmod('PF_implementacion_modeloL2',ODE_ss,u);

SS_desv = (F_0-ODE_ss)';

L_desv = (L_out-L_out(1))';

%% Análisis de respuestas
Sys_ss = ss(A1,B1,C1,D1);   %Objeto-espacio de estados % D,'StateName',{'Ca';'Cb'},'InputName',{'D';'Caf'},'OutputName',{'Ca';'Cb'}
Sys_tf = tf(Sys_ss);        %Transformación a función de tranferencia

%% Reducción

Sys_min = minreal(Sys_ss);

mSys_min = sminreal(Sys_ss);

rSys_min = balred(Sys_tf,2);



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

t_sim = 9000;

t_sp_change = 5000;

%% Flujo

%% SP y DR

L_ss_ini = L_ss*0.5;

L_ss_fin = L_ss;

%%
% out = sim('PF_implementacion_PID_Flujo');

%save('out_PID_C4.mat','out')


% F_in_t = out.Inputs.Data(:,1);
% S1_in_t = out.Inputs.Data(:,3);
% 
% L_ODE = out.Outputs_Reactor.Data(:,5);
% 
% L_error = out.L_set.Data(:,1);
% 
% time_inputs = out.Inputs.Time;
% 
% SP_changes = out.SP_track.Data(:,1);


% figure(1)
% subplot(2,1,1)
% plot(time_inputs,L_ODE,'LineWidth',2.5)
% hold on
% plot(time_inputs,SP_changes,'r--','LineWidth',1.5)
% 
% ylabel('L')
% title('Seguimiento de Set Point')
% legend('L(t)','L_{set}','Location','northwest')
% grid on
% 
% subplot(2,1,2)
% 
% yyaxis left
% plot(time_inputs,F_in_t,'LineWidth',2)
% ylabel('F_{in}')
% 
% yyaxis right
% stairs(time_inputs,S1_in_t,'LineWidth',2)
% ylabel('S_{1,in} (g/L)')
% 
% xlabel('Tiempo')
% title('Perturbación aplicada')
% grid on
% 
% legend('F_{in}','S_{1,in}','Location','best')
% 
% set(gcf,'Color','w')

% exportgraphics(gcf,...
%     'PID_STP_DF_Flujo_C.png',...
%     'Resolution',600)

% v = VideoWriter('PID_STP_DF_Flujo_C.mp4','MPEG-4');
% v.FrameRate = 20;
% open(v)
% 
% fig = figure('Position',[100 100 1000 700]);
% 
% for k = 1:length(time_inputs)
% 
%     clf
% 
%     %% Subplot superior
%     subplot(2,1,1)
% 
%     plot(time_inputs(1:k),L_ODE(1:k),...
%         'LineWidth',2.5)
%     hold on
% 
%     plot(time_inputs(1:k),SP_changes(1:k),...
%         'r--','LineWidth',1.5)
% 
%     ylabel('L')
%     title('Seguimiento de Set Point')
% 
%     legend('L(t)','L_{set}',...
%         'Location','northwest')
% 
%     grid on
%     xlim([min(time_inputs) max(time_inputs)])
% 
%     ymin = min([L_ODE(:); SP_changes(:)]);
%     ymax = max([L_ODE(:); SP_changes(:)]);
% 
%     ylim([0.95*ymin 1.05*ymax])
% 
%     %% Subplot inferior
%     subplot(2,1,2)
% 
%     yyaxis left
% 
%     plot(time_inputs(1:k),...
%          F_in_t(1:k),...
%          'LineWidth',2.5)
% 
%     ylabel('F_{in}')
% 
%     yminL = min(F_in_t);
%     ymaxL = max(F_in_t);
% 
%     if yminL == ymaxL
%         ylim([0.9*yminL 1.1*ymaxL])
%     else
%         ylim([0.95*yminL 1.05*ymaxL])
%     end
% 
%     yyaxis right
% 
%     stairs(time_inputs(1:k),...
%            S1_in_t(1:k),...
%            'LineWidth',2.5)
% 
%     ylabel('S_{1,in} (g/L)')
% 
%     yminR = min(S1_in_t);
%     ymaxR = max(S1_in_t);
% 
%     if yminR == ymaxR
%         ylim([0.9*yminR 1.1*ymaxR])
%     else
%         ylim([0.95*yminR 1.05*ymaxR])
%     end
% 
%     xlabel('Tiempo')
%     title('Perturbación aplicada')
% 
%     grid on
%     xlim([min(time_inputs) max(time_inputs)])
% 
%     legend({'F_{in}','S_{1,in}'},...
%            'Location','northwest')
% 
%     ax = gca;
%     ax.YAxis(1).Color = 'b';
%     ax.YAxis(2).Color = 'r';
% 
%     %% Guardar frame
%     frame = getframe(fig);
%     writeVideo(v,frame);
% 
% end
% 
% close(v)


%%
% out1 = sim('PF_implementacion_PID_Intensidad');
% 
% I0_in_t = out1.Inputs.Data(:,7);
% S1_in_t = out1.Inputs.Data(:,3);
% 
% L_ODE = out1.Outputs_Reactor.Data(:,5);
% 
% L_error = out1.L_set.Data(:,1);
% 
% time_inputs = out1.Inputs.Time;
% 
% SP_changes = out1.SP_track.Data(:,1);
% 
% figure(1)
% subplot(2,1,1)
% plot(time_inputs,L_ODE,'LineWidth',2.5)
% hold on
% plot(time_inputs,SP_changes,'r--','LineWidth',1.5)
% 
% ylabel('L')
% title('Seguimiento de Set Point')
% legend('L(t)','L_{set}','Location','northwest')
% grid on
% 
% subplot(2,1,2)
% 
% yyaxis left
% plot(time_inputs,I0_in_t,'LineWidth',2)
% ylabel('I_0')
% 
% yyaxis right
% stairs(time_inputs,S1_in_t,'LineWidth',2)
% ylabel('S_{1,in} (g/L)')
% 
% xlabel('Tiempo')
% title('Perturbación aplicada')
% grid on
% 
% legend('I_0','S_{1,in}','Location','best')
% 
% set(gcf,'Color','w')
% 
% exportgraphics(gcf,...
%     'PID_STP_DF_Intensidad_C.png',...
%     'Resolution',600)


%%
% out2 = sim('PF_implementacion_PID_S1');
% 
% 
% S1_in_t = out2.Inputs.Data(:,3);
% I0_in_t = out2.Inputs.Data(:,7);
% 
% L_ODE = out2.Outputs_Reactor.Data(:,5);
% 
% L_error = out2.L_set.Data(:,1);
% 
% time_inputs = out2.Inputs.Time;
% 
% SP_changes = out2.SP_track.Data(:,1);
% 
% figure(1)
% subplot(2,1,1)
% plot(time_inputs,L_ODE,'LineWidth',2.5)
% hold on
% plot(time_inputs,SP_changes,'r--','LineWidth',1.5)
% 
% ylabel('L')
% title('Seguimiento de Set Point')
% legend('L(t)','L_{set}','Location','northwest')
% grid on
% 
% subplot(2,1,2)
% 
% yyaxis left
% plot(time_inputs,S1_in_t,'LineWidth',2)
% ylabel('S_{1,in} (g/L)')
% 
% yyaxis right
% stairs(time_inputs,I0_in_t,'LineWidth',2)
% ylabel('I_0')
% 
% xlabel('Tiempo')
% title('Perturbación aplicada')
% grid on
% 
% legend('S_{1,in}','I_0','Location','best')
% 
% set(gcf,'Color','w')

% exportgraphics(gcf,...
%     'PID_STP_DF_S1_C.png',...
%     'Resolution',600)

%%





%% Ruido sensor

% out_F_R = sim('PF_implementacion_PID_Flujo_Ruido');
% 
% F_in_t_R = out_F_R.Inputs.Data(:,1);
% S1_in_t_R = out_F_R.Inputs.Data(:,3);
% 
% L_ODE_R = out_F_R.Outputs_Reactor.Data(:,5);
% 
% L_error_R = out_F_R.L_set.Data(:,1);
% 
% time_inputs_R = out_F_R.Inputs.Time;
% 
% SP_changes_R = out_F_R.SP_track.Data(:,1);
% 
% figure(1)
% subplot(2,1,1)
% 
% plot(time_inputs_R,L_ODE_R,'LineWidth',2.5)
% hold on
% plot(time_inputs_R,SP_changes_R,'r--','LineWidth',1.5)
% 
% ylabel('L')
% title('Seguimiento de Set Point')
% legend('L(t)','L_{set}','Location','northwest')
% grid on
% 
% subplot(2,1,2)
% 
% stairs(time_inputs_R,S1_in_t_R,'LineWidth',2)
% 
% ylabel('S_{1,in} (g/L)')
% xlabel('Tiempo')
% title('Perturbación aplicada')
% grid on
% 
% set(gcf,'Color','w')
% 
% exportgraphics(gcf,...
%     'PID_STP_DF_Flujo_Ruido.png',...
%     'Resolution',600)




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
