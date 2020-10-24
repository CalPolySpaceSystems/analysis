"""This module estimates altitude for a rocket with constant thrust and constant drag coefficient"""

from math import pi
from scipy.integrate import solve_ivp
from ambiance import Atmosphere
import matplotlib.pyplot as plt

THRUST = 2560 # N
IMPULSE = 40960 # N*s
PROP_MASS = 17.1 # kg
DRY_MASS = 95 * .452592 # kg
N2_MASS = 3 # kg
DIAMETER = .22 # m
DRAG_COEFFICIENT = .5 # unitless
LAUNCH_ELEVATION = 1400 # m, Spaceport America

A_C = pi/4 * DIAMETER ** 2
WET_MASS = DRY_MASS + PROP_MASS


def thrust_curve(time):
    """Returns a flat thrust curve."""
    burn_time = IMPULSE/THRUST
    if time < burn_time:
        current_thrust = THRUST
    else:
        current_thrust = 0.0
    return current_thrust

def mass_curve(time):
    """Returns a mass curve for a flat thrust curve"""
    burn_time = IMPULSE/THRUST
    if time < burn_time:
        mass = DRY_MASS + (1 - time/burn_time) * PROP_MASS
    else:
        mass = DRY_MASS - N2_MASS
    return mass

def drag(height, velocity):
    """Returns drag force from altitude and velocity"""
    atm = Atmosphere(height)
    return DRAG_COEFFICIENT * A_C * .5 * atm.density * velocity * abs(velocity)

def flight(t, y):
    """ODE defining the trajectory of the rocket"""
    mass = mass_curve(t)
    d_t = [y[1], thrust_curve(t) / mass - 9.807 - drag(y[0], y[1]) / mass]
    return d_t

solution = solve_ivp(flight, [0, 90], [LAUNCH_ELEVATION, 0], rtol=1e-10, atol=1e-12)
plt.plot(solution.t, solution.y[0]-LAUNCH_ELEVATION)
plt.xlabel('Time (s)')
plt.ylabel('Altitude (m)')
plt.show()
