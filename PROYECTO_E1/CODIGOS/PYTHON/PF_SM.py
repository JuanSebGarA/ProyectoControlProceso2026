#A simulation model for bio - manufacturing process development and analysis
# Description :
# This code contains a model that solves the dynamical and steady state of
# CSTR bioreactor for microalgae and a tubular PBR . Assumptions and
# derivation of the model equations are shown in the Thesis report .
# The model considers constant volume for the reactors
import numpy as np
from scipy.integrate import solve_ivp # To solve the initial value problem
import matplotlib.pyplot as plt

# Parameters from the model equations
mu_m = 0.0418 # Cell specific growth rate [1/ h ]
q_0 = 0.0196 # Minimum nitrogen quota [ g / g ]
l_0 = 0.0006 # Minimum lipid quota [ g / g ]
K_S1 = 0.1002 # Half saturation constant of C substrate for growth [ g / L ]
K_I = 66.5337 # Half saturation constant for light growth [ umol / m2 / s ]
K_I2 = 100.0 # Inhibition constant for light growht [ umol / m2 / s ]
m_s = 0.0167 # Maintenance coefficent [1/ h ]
rho_m = 0.1197 # Maximum nitrogen uptake rate into the cell [1/ h ]
K_N = 0.5793 # Half - saturation constant of glycine for nitrogen uptake [ g / L ]
q_m = 0.2109 # Maximum nitrogen quota [ g / g ]
pi_m = 0.0762 # Maximum lipid production rate [1/ h ]
K_L = 12.5596 # Half - saturation constant of C substrate for lipid production [ g / L ]
l_m = 0.6995 # Maximum lipid quota [ g / g ]
Y_xs = 0.9597 # Substrate - to - biomass yield coefficient [ g / g ]
Y_ls = 0.1908 # Substrate - to - lipid yield coefficient [ g / g ]
alpha = 25.4804
B = 0.04

model_params = { 'mu_m': mu_m , 'q_0': q_0 , 'l_0': l_0 , 'K_S1': K_S1 , 
                'K_I': K_I , 'K_I2': K_I2 , 'm_s': m_s , 'rho_m': rho_m , 
                'K_N': K_N , 'q_m': q_m , 'pi_m' : pi_m , 'K_L': K_L , 'l_m': l_m , 
                'Y_xs': Y_xs , 'Y_ls': Y_ls , 'alpha': alpha , 'B': B }

# KINETIC EXPRESSIONS : RATE EQUATIONS

def f_kinetic_biomass (q , l , S_1 , I , model_params ):
	mu_m , q_0 , l_0 , K_S1 , K_I , K_I2 = [ model_params[p] for p in [ 'mu_m', 'q_0' , 'l_0' , 'K_S1', 'K_I' , 'K_I2' ]]
	mu = mu_m * (1 - q_0 / q ) * (1 - l_0 / l ) * ( S_1 /( K_S1 + S_1 ) ) * ( I /( K_I + I +( I**2) / K_I2 ) )
	mu = mu *( mu >=0) *( q >=0) *( l >=0) *( q >= q_0 ) *( l >= l_0 ) *( S_1 >=0) 
	# Zero if exceeding max or min limits
	# if mu > mu_m
	# mu = mu_m
	return mu


# Nitrogen uptake rate into the cell :
def f_kinetic_nitrogen_uptake (q , S_2 , model_params ):
    rho_m , q_m , q_0 , K_N = [ model_params[p] for p in [ 'rho_m' , 'q_m' , 'q_0' , 'K_N' ]]
    rho = rho_m * (( q_m - q ) /( q_m - q_0 ) ) * ( S_2 /( K_N + S_2 ) )
    rho = rho *( rho >=0) *( S_2 >=0) *( q >=0) *( S_2 >=0) # Zero if exceeding max or min limits
	# if rho > rho_m :
	# rho = rho_m
    return rho

# Lipid production rate :
def f_kinetic_lipid_production ( S_1 , q , l , model_params ) :
	pi_m , K_L , l_m , q_m = [ model_params[p] for p in [ 'pi_m' , 'K_L' , 'l_m' , 'q_m' ]]
	pi = pi_m * ( S_1 /( K_L + S_1 ) ) * (1 - q ) * (( l_m - l ) / l_m )
	pi = pi *( pi >=0) *( S_1 >=0) *( q >=0) *( l >=0) *( q < q_m ) *( l < l_m ) # Zero if exceeding max or min limits
	# if pi > pi_m :
	# pi = pi_m
	return pi

# Average light irradiance in the reactor :
def f_culture_irradience ( XX , I0 , z , model_params ) :
	alpha , B = [ model_params [ p ] for p in [ 'alpha' , 'B' ]]
	BX = B * XX # Chlorophyll concentration [ g / L ]
	#I = 2* I0 /( 1000*alpha * BX * z **2) *( z -1/( 1000*alpha * BX ) + np . exp ( - 1000*alpha * BX * z ) /( 1000*alpha * BX ) )
	I = 2* I0 /( alpha * BX * z **2) *( z -1/( alpha * BX ) + np . exp ( -alpha * BX * z ) /( alpha * BX ) )
	return I

# FUNCTIONS OF THE MODELS
def f_cstr_model_time (t , x , u , par , model_params ):
	
	# Initial concentrations of the species
	
	X = x [0]
	S_1 = x [1]
	S_2 = x [2]
	Q = x [3]
	L = x [4]
	
	# Calculation of total biomass and quotas
	
	XX = X + Q + L # Total biomass concentration g / L
	l = L / XX # Lipid quota [ g / g ]
	q = Q / XX # Nitrogen quota [ g / g ]
	
	# Inputs
	
	F_in = u [0] # L / h
	I0 = u [1] # umol / m2 / s
	X_in = u [2]
	S_1_in = u [3]
	S_2_in = u [4]
	Q_in = u [5]
	L_in = u [6]
	
	# Reactor geometry parameters
	
	VOL = par [0]
	R = par [1]
	
	# Overflow control
	
	F_out = F_in
	
	# Kinetic expressions
	
	I = f_culture_irradience ( XX , I0 , R , model_params )
	mu = f_kinetic_biomass (q , l , S_1 , I , model_params )
	rho = f_kinetic_nitrogen_uptake (q , S_2 , model_params )
	pi = f_kinetic_lipid_production ( S_1 , q , l , model_params )
	
	Y_xs , Y_ls , m_s = [ model_params [ p ] for p in [ 'Y_xs' , 'Y_ls' , 'm_s']]
	
	# Differential equations of the model
	
	dXdt = ( mu * X ) + X_in * F_in / VOL - X * F_out / VOL
	dS_1dt = ( -1/ Y_xs * mu *X -1/ Y_ls * pi * X ) - m_s * X + S_1_in * F_in / VOL - S_1 * F_out / VOL
	dS_2dt = ( - rho * X ) + S_2_in * F_in / VOL - S_2 * F_out / VOL
	dQdt = ( rho *X - mu * Q ) + Q_in * F_in / VOL - Q * F_out / VOL
	dLdt = ( pi *X - mu * L ) + L_in * F_in / VOL - L * F_out / VOL
	dVdt = F_in - F_out
	
	# Storing the differential terms in a vector
	
	xdot = np.zeros (5)
	xdot [0] = dXdt
	xdot [1] = dS_1dt
	xdot [2] = dS_2dt
	xdot [3] = dQdt
	xdot [4] = dLdt
	# xdot [5] = dVdt
	
	return xdot



Vol = 270
R = 0.15 # 1.5 en el caso de la mejor producción

par = [Vol, R]

F_in   = 1.5   # Caudal de entrada      [L/h]
I0     = 200.0  # Irradiancia incidente  [umol/m2/s]
X_in   = 0.0    # Biomasa entrante       [g/L]
S1_in  = 60.0    # Sustrato C entrante    [g/L]
S2_in  = 20.0    # Nitrógeno entrante     [g/L]
Q_in   = 0.0   # Cuota N entrante       [g/L]
L_in   = 0.0    # Lípidos entrantes      [g/L]

u = [F_in, I0, X_in, S1_in, S2_in, Q_in, L_in]

x0 = [0.25,   # X  - Biomasa estructural   [g/L]
      0,   # S1 - Sustrato carbono      [g/L]
      0,   # S2 - Nitrógeno             [g/L]
      0.025,  # Q  - Cuota nitrógeno       [g/L]
      0.025]  # L  - Lípidos               [g/L]

# ── Horizonte de simulación ─────────────────────────────────────────────────
t_start = 0.0
t_end   = 1250.0          # horas
t_eval  = np.linspace(t_start, t_end, 1000)

# ── Resolución de la ODE ────────────────────────────────────────────────────
sol = solve_ivp(
    fun      = lambda t, x: f_cstr_model_time(t, x, u, par, model_params),
    t_span   = (t_start, t_end),
    y0       = x0,
    method   = 'RK45',       # Runge-Kutta 4(5) — bueno para sistemas suaves
    t_eval   = t_eval,
    rtol     = 1e-6,
    atol     = 1e-8
)

t     = sol.t
x_out = sol.y[0]   # X  - Biomasa estructural  [g/L]
Q_out = sol.y[3]   # Q  - Cuota de nitrógeno   [g/L]
L_out = sol.y[4]   # L  - Lípidos              [g/L]
X_out = x_out + Q_out + L_out   # Biomasa total [g/L]

# ── Gráfico (replica exacta del figure(1) de MATLAB) ───────────────────────
plt.figure(1, figsize=(9, 5))

plt.plot(t, X_out, linewidth=2, color='black', label='$X_{total}$')
plt.plot(t, Q_out, linewidth=2, color='red',   label='$Q$')
plt.plot(t, L_out, linewidth=2, color='green', label='$L$', linestyle = 'dashed')
plt.plot(t, x_out, linewidth=2, color='blue',  label='$X$')

plt.xlabel('Tiempo [h]', fontsize=12)
plt.ylabel('Concentración [g/L]', fontsize=12)
plt.title(f'Modelo Dinamico – $F_{{in}}$ = {F_in} L/h', fontsize=13)
plt.legend(fontsize=11)
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.show()
