function [sys,x0,str,ts,simStateCompliance] = PF_Sfunction_L(t,x,u,flag)
%SFUNTMPL General MATLAB S-Function Template
%   With MATLAB S-functions, you can define you own ordinary differential
%   equations (ODEs), discrete system equations, and/or just about
%   any type of algorithm to be used within a Simulink block diagram.
%
%   The general form of an MATLAB S-function syntax is:
%       [SYS,X0,STR,TS,SIMSTATECOMPLIANCE] = SFUNC(T,X,U,FLAG,P1,...,Pn)
%
%   What is returned by SFUNC at a given point in time, T, depends on the
%   value of the FLAG, the current state vector, X, and the current
%   input vector, U.
%
%   FLAG    RESULT              DESCRIPTION
%   -----   ------              --------------------------------------------
%   0       [SIZES,X0,STR,TS]   Initialization, return system sizes in SYS,
%                               initial state in X0, state ordering strings
%                               in STR, and sample times in TS.
%   1       DX                  Return continuous state derivatives in SYS.
%   2       DS                  Update discrete states SYS = X(n+1)
%   3       Y                   Return outputs in SYS.
%   4       TNEXT               Return next time hit for variable step sample
%                               time in SYS.
%   5                           Reserved for future (root finding).
%   9       []                  Termination, perform any cleanup SYS=[].
%
%
%   The state vectors, X and X0 consists of continuous states followed
%   by discrete states.
%
%   Optional parameters, P1,...,Pn can be provided to the S-function and
%   used during any FLAG operation.
%
%   When SFUNC is called with FLAG = 0, the following information
%   should be returned:
%
%     SYS(1) = Number of continuous states.
%     SYS(2) = Number of discrete states.
%     SYS(3) = Number of outputs.
%     SYS(4) = Number of inputs.
%              Any of the first four elements in SYS can be specified
%              as -1 indicating that they are dynamically sized. The
%              actual length for all other flags will be equal to the
%              length of the input, U.
%              SYS(5) = Reserved for root finding. Must be zero.
%              SYS(6) = Direct feedthrough flag (1=yes, 0=no). The s-function
%              has direct feedthrough if U is used during the FLAG=3
%              call. Setting this to 0 is akin to making a promise that
%              U will not be used during FLAG=3. If you break the promise
%              then unpredictable results will occur.
%     SYS(7) = Number of sample times. This is the number of rows in TS.
%
%
%     X0     = Initial state conditions or [] if no states.
%
%     STR    = State ordering strings which is generally specified as [].
%
%     TS     = An m-by-2 matrix containing the sample time
%               (period, offset) information. Where m = number of sample
%               times. The ordering of the sample times must be:
%
%               TS = [0         0,      : Continuous sample time.
%                     0         1,      : Continuous, but fixed in minor step
%                                         sample time.
%                     PERIOD    OFFSET, : Discrete sample time where
%                                         PERIOD > 0 & OFFSET < PERIOD.
%                      -2       0];     : Variable step discrete sample time
%                                         where FLAG=4 is used to get time of
%                                         next hit.
%
%               There can be more than one sample time providing
%               they are ordered such that they are monotonically
%               increasing. Only the needed sample times should be
%               specified in TS. When specifying more than one
%               sample time, you must check for sample hits explicitly by
%               seeing if
%                  abs(round((T-OFFSET)/PERIOD) - (T-OFFSET)/PERIOD)
%               is within a specified tolerance, generally 1e-8. This
%               tolerance is dependent upon your model's sampling times
%               and simulation time.
%
%               You can also specify that the sample time of the S-function
%               is inherited from the driving block. For functions which
%               change during minor steps, this is done by
%               specifying SYS(7) = 1 and TS = [-1 0]. For functions which
%               are held during minor steps, this is done by specifying
%               SYS(7) = 1 and TS = [-1 1].
%
%       SIMSTATECOMPLIANCE = Specifices how to handle this block when saving and
%                            restoring the complete simulation state of the
%                            model. The allowed values are: 'DefaultSimState',
%                            'HasNoSimState' or 'DisallowSimState'. If this value
%                            is not speficified, then the block's compliance with
%                            simState feature is set to 'UknownSimState'.
%
%       Copyright 1990-2010 The MathWorks, Inc.
%
% The following outlines the general structure of an S-function.

switch flag,
    %%%%%%%%%%%%%%%%%%
    % Initialization %
    %%%%%%%%%%%%%%%%%%
    case 0,
        [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes;
        %%%%%%%%%%%%%%%
        % Derivatives %
        %%%%%%%%%%%%%%%
    case 1,
        sys=mdlDerivatives(t,x,u);
        %%%%%%%%%%
        % Update %
        %%%%%%%%%%
    case 2,
        sys=mdlUpdate(t,x,u);
        %%%%%%%%%%%
        % Outputs %
        %%%%%%%%%%%
    case 3,
        sys=mdlOutputs(t,x,u);
        %%%%%%%%%%%%%%%%%%%%%%%
        % GetTimeOfNextVarHit %
        %%%%%%%%%%%%%%%%%%%%%%%
    case 4,
        sys=mdlGetTimeOfNextVarHit(t,x,u);
        %%%%%%%%%%%%%
        % Terminate %
        %%%%%%%%%%%%%
    case 9,
        sys=mdlTerminate(t,x,u);
        %%%%%%%%%%%%%%%%%%%%
        % Unexpected flags %
        %%%%%%%%%%%%%%%%%%%%
    otherwise
        DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));
end
% end sfuntmpl
%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes
%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded. This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;
sizes.NumContStates  = 5; % Número de estados continuos
sizes.NumDiscStates  = 0; % Número de estados discretos
sizes.NumOutputs     = 5; % Número de salidas
sizes.NumInputs      = 6; % Número de entradas
sizes.DirFeedthrough = 1; % Direct feedthrough flag (1=yes, 0=no)
                          % Para indicar si las entradas afectan la salida
sizes.NumSampleTimes = 1; % Número de filas en el vector de tiempo de muestreo
                          % (por lo menos debe ser 1)
sys = simsizes(sizes);

%
% initialize the initial conditions

x0 = [14.7023102022666;	0.229849249096564;	11.7210309830145;	4.04050938037401;	1.09061064635207];
%x0 = [0.25; 0; 0; 0.025; 0.025];
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts = [0 0]; % Sistema continuo

%   Specify the block simStateCompliance. The allowed values are:
%       'UnknownSimState', < The default setting; warn and assume DefaultSimState
%       'DefaultSimState', < Same sim state as a built-in block
%       'HasNoSimState',   < No sim state
%       'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes
%

%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u)
%Entradas
% Parametros Entrada
I0 = 200;
R = 2.0;
alpha = 0.254804;
B = 0.04;
V = 270;

F_in = u(1); x_in = u(2); S1_in = u(3); S2_in = u(4);
Q_in = u(5); L_in = u(6); 

% Parámetros maximos
max = [0.0418 0.2109 0.6995 0.1197 0.0762]; % mu,q,l,rho,pi
% Parámetros minimos
min = [0.0196 0.0006]; % q,l
% Parámetros modelo
param_mod = [0.9597 0.1908 0.0167 0.1002 0.579 12.5596 66.5337 100]; % Yxs, Yls, ms, mu_max, qmin, qmax, lmin, lmax, rho_max, pi_max, KS1, KN, KL, KI, KI2
%Salidas
% x_f S1 S2 Q L
x_f = x(1);       S1 = x(2);
S2 = x(3);      Q = x(4);
L = x(5);

%Parametros Max-Min
mu_max = max(1); q_max = max(2); l_max = max(3); rho_max = max(4); pi_max = max(5);
q_min = min(1); l_min = min(2);


% Parametros modelos
Yxs = param_mod(1); Yls = param_mod(2); ms = param_mod(3); K_S1 = param_mod(4);
K_N = param_mod(5); K_L = param_mod(6); K_I = param_mod(7); K_I2 = param_mod(8);

% Totales
X = x_f + Q + L;
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
dx_dt = mu * x_f + (F_in/V)*(x_in - x_f);
dS1_dt = -(1/Yxs)*(mu * x_f) - (1/Yls)*(pi*x_f) - ms*x_f + (F_in/V)*(S1_in-S1);
dS2_dt = - rho*x_f + (F_in/V)*(S2_in-S2);
dQ_dt = rho*x_f - mu * Q + (F_in/V) * (Q_in - Q);
dL_dt = pi * x_f - mu * L + (F_in/V) * (L_in - L);


sys = [dx_dt; dS1_dt; dS2_dt; dQ_dt; dL_dt];
% end mdlDerivatives
%

%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%

function sys=mdlUpdate(t,x,u)
sys = []; % Ecuaciones para los estados discretos (no se usan)
% end mdlUpdate
%

%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)
mc = eye(length(x)); %matriz C del espacio de estados
sys = mc*x;        % Salidas del sistema continuo (multiplicación matricial)

% end mdlOutputs
%

%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block. Note that the result is
% absolute time. Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)
sampleTime = 1; % Example, set the next hit to be one second later.
sys = t + sampleTime;
% end mdlGetTimeOfNextVarHit
%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)
sys = [];
% end mdlTerminate
%=============================================================================