"""This module estimates altitude for a rocket with constant thrust and constant drag coefficient"""

from math import pi
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from ambiance import Atmosphere
import matplotlib.pyplot as plt
from motor_performance.rpa_wrapper import RPASim
from motor_performance.transient_analysis import RocketPerfomance

## Prop system conditions
BURN_TIME = 15.0 # s, burn time
IMPULSE = 40960 # N*s, upper limit on impulse
M_DOT_OX = 1.02 # kg/s, oxidizer mdot
ID = 2.15 # in, initial port diamater
OD = 3.4 # in, outer grain diameter
L = 21.607 # in, grain length
D_THROAT = 1.0 * .0254 # m, throat diameter
EXPANSION_RATIO = 6.25
EXIT_AREA = pi/4 * D_THROAT**2 * EXPANSION_RATIO # m^2

# fuel characteristics
RHO_FUEL = 964          # fuel density, kg/m3
A = .417
N = .347      # ballistic coefficients, G_ox = g/(cm2*s), rdot = mm/s
# For more info on what these mean, see Sutton

## Other rocket conditions
DRY_MASS = 95 * .452592 # kg, mass without N2 or propellant
N2_MASS = 3 # kg
DIAMETER = .22 # m, rocket diameter
DRAG_COEFFICIENT = .5 # unitless
LAUNCH_ELEVATION = 1400 # m, Spaceport America

# Derived paramaters
PROP_MASS = BURN_TIME * M_DOT_OX + pi/4 * ((OD*.0254)**2 - (ID*.0254)**2) * L*.0254 * RHO_FUEL
A_C = pi/4 * DIAMETER ** 2
WET_MASS = DRY_MASS + PROP_MASS

C_D_CURVE = np.loadtxt(fname='altitude_approximation/Cd vs Flight Speed.txt', skiprows=1)
c_d_func = interp1d(C_D_CURVE[:,0], C_D_CURVE[:,2], fill_value='extrapolate')
plt.figure()
plt.plot(C_D_CURVE[:,0], c_d_func(C_D_CURVE[:,0]))
plt.plot(C_D_CURVE[:,0], C_D_CURVE[:,2])
plt.show()

# Solve for grain regresseion characteristics
regression_rate = lambda t, y: .1 * A * (M_DOT_OX*1000 / (pi * y**2))**N
sol = solve_ivp(regression_rate, [0, BURN_TIME], [ID * 1.27], rtol=1e-10, atol=1e-12)
r_dot = 10 * (sol.y[0][1:] - sol.y[0][:-1]) / (sol.t[1:] - sol.t[:-1]) # mm/s, regression rate
m_dot_fuel = 2 * pi**(1-N) * RHO_FUEL * L * .0254 * A * (M_DOT_OX * 1000)**N * sol.y[0]**(1 - 2*N) / 100000 # kg/s, vector over time
o_f_ratio = M_DOT_OX / m_dot_fuel # vector of O/F over time

hm9 = RPASim(path='../rpa_lib')
hm9.sim_config(p_c=40.0, fuel='HTPB', oxidizer='N2O(L),298.15K', o_f_ratio=o_f_ratio[0])
hm9.model_chamber(d_throat=[D_THROAT, 'm'])
hm9.model_nozzle_flow(exit_condition='area ratio', exit_condition_value=EXPANSION_RATIO)

performance = RocketPerfomance(rocket=hm9, m_dot_fuel=m_dot_fuel, m_dot_ox=M_DOT_OX)

expansion_ratios = np.linspace(6, 8, 11)
exit_areas = pi/4 * D_THROAT**2 * expansion_ratios
v_exit_poly = np.polyfit(x=sol.t, y=performance.v_exit_curve, deg=5)
m_dot_poly = np.polyfit(x=sol.t, y=m_dot_fuel+M_DOT_OX, deg=5)
p_exit_poly = np.polyfit(x=sol.t, y=performance.p_exit_curve, deg=5)

def get_thrust_at_altitude(time, p_ambient, exit_area=[EXIT_AREA]):
    """Returns a flat thrust curve."""
    v_exit = np.polyval(v_exit_poly, time)
    p_exit = np.polyval(p_exit_poly, time)
    m_dot = np.polyval(m_dot_poly, time)
    if time < BURN_TIME:
        return m_dot * v_exit + (p_exit * 10**5 - p_ambient) * exit_area
    else:
        return 0

def drag(velocity, density):
    """Returns drag force from altitude and velocity"""
    return c_d_func(abs(velocity)) * A_C * .5 * density * velocity * abs(velocity)

def flight(t, y):
    """ODE defining the trajectory of the rocket"""
    d_t = np.zeros(4)
    atm = Atmosphere(y[0])
    d_t[0] = y[1]
    d_t[3] = get_thrust_at_altitude(t, atm.pressure, [exit_area]) if y[3] < IMPULSE else 0
    d_t[1] = d_t[3] / y[2] - 9.807 - drag(y[1], atm.density) / y[2]
    d_t[2] = -np.polyval(m_dot_poly, t) if d_t[3] > 0 else 0.0
    return d_t

solutions = []
for i in range(len(exit_areas)):
    exit_area = exit_areas[i]
    performance.generate_exit_curves(exit_area=exit_area)
    v_exit_poly = np.polyfit(x=sol.t, y=performance.v_exit_curve, deg=5)
    p_exit_poly = np.polyfit(x=sol.t, y=performance.p_exit_curve, deg=5)
    solutions.append(solve_ivp(flight, [0, 50], [LAUNCH_ELEVATION, 0, WET_MASS, 0], rtol=1e-10, atol=1e-12))

legend_labels = []
for ratio in expansion_ratios:
    legend_labels.append(f'AR = {ratio}')


plt.figure()
for solution in solutions:
    plt.plot(solution.t, solution.y[0]-LAUNCH_ELEVATION)
plt.xlabel('Time (s)')
plt.ylabel('Altitude (m)')
plt.legend(legend_labels)

plt.figure()
for solution in solutions:
    plt.plot(solution.t, solution.y[1])
plt.xlabel('Time (s)')
plt.ylabel('Velocity (m/s)')
plt.legend(legend_labels)

plt.figure()
for solution in solutions:
    plt.plot(solution.t, solution.y[2])
plt.xlabel('Time (s)')
plt.ylabel('Mass (kg)')
plt.legend(legend_labels)

plt.figure()
for solution in solutions:
    plt.plot(solution.t, solution.y[3])
plt.xlabel('Time (s)')
plt.ylabel('Impulse (N*s)')
plt.legend(legend_labels)

plt.show()
