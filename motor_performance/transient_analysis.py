from math import pi
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
# from rpa_wrapper import HM9
#
# # Designed system to test
# M_DOT_OX = 1.02 # kg/s, oxidizer mdot
# ID = 2.15 # in, initial port diamater
# L = 21.607 # in, grain length
# D_THROAT = 1.064 # in, throat diameter
# EXPANSION_RATIO = 6.25
# BURN_TIME = 14 # s, burn time
#
# # fuel characteristics
# RHO_FUEL = 964          # fuel density [kg/m3]
# A = .417; N = .347      # ballistic coefficients [G_ox = g/(cm2*s), rdot = mm/s]
# # For more info on what these mean, see Sutton
#
# regression_rate = lambda t, y: .1 * A * (M_DOT_OX*1000 / (pi * y**2))**N
# sol = solve_ivp(regression_rate, [0, BURN_TIME], [ID * 1.27], rtol=1e-10, atol=1e-12)
# r_dot = 10 * (sol.y[0][1:] - sol.y[0][:-1]) / (sol.t[1:] - sol.t[:-1]) # mm/s, regression rate
#
# m_dot_fuel = 2 * pi**(1-N) * RHO_FUEL * L * .0254 * A * (M_DOT_OX * 1000)**N * sol.y[0]**(1 - 2*N) / 100000 # kg/s
# o_f_ratio = M_DOT_OX / m_dot_fuel
#
# hybrid = HM9(m_dot_fuel[0], M_DOT_OX)

class RocketPerfomance():
    def __init__(self, rocket, m_dot_fuel, m_dot_ox, initial_pc=35.0):
        print('init RocketPerfomance')
        self.p_c_curve = np.ones(len(m_dot_fuel))
        self.thrust_curve = np.ones(len(m_dot_fuel))
        self.isp_curve = np.ones(len(m_dot_fuel))
        self.v_exit_curve = np.ones(len(m_dot_fuel))
        self.p_exit_curve = np.ones(len(m_dot_fuel))
        self.throat_velocity_curve = np.ones(len(m_dot_fuel))
        self.rocket = rocket

        for i in range(len(m_dot_fuel)):
            def _chamber_pressure_function(guess):
                m_dot_f = m_dot_fuel[i]
                self.rocket.rebuild(p_c=guess, m_dot_fuel=m_dot_f, m_dot_ox=m_dot_ox)
                return self.rocket.m_dot - m_dot_f - m_dot_ox

            if self.p_c_curve.all() == np.empty(len(m_dot_fuel)).all():
                guess = initial_pc
            else:
                guess = self.p_c_curve[-1]

            p_c = fsolve(_chamber_pressure_function, guess)
            self.p_c_curve[i] = p_c
            self.thrust_curve[i] = self.rocket.thrust
            self.isp_curve[i] = self.rocket.isp
            self.throat_velocity_curve[i] = self.rocket.throat_velocity

    def generate_exit_curves(self, exit_area):
        k = 1.2
        for i in range(len(self.p_c_curve)):
            p_c = self.p_c_curve[i]
            def _exit_pressure_function(guess):
                # eq 3-25 in Sutton
                return ((k+1)/2) ** (1/(k-1)) * (guess/p_c) ** (1/k) * ((k+1)/(k-1) * (1 - (guess/p_c) ** ((k-1)/k))) ** (.5) * exit_area - pi/4 * self.rocket.d_throat**2
            p_e = fsolve(_exit_pressure_function, 1.0)
            print(p_e)
            self.p_exit_curve[i] = p_e
            v_e = self.throat_velocity_curve[i] * ((k+1)/(k-1) * (1 - (p_e/p_c) ** ((k-1)/k))) ** (.5) # m/s
            print(v_e)
            self.v_exit_curve[i] = v_e
