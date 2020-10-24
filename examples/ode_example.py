"""An example of how to solve a system of ODEs. Uses the lorenz equation"""

import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

SIGMA = 10
R = 28
B = 8/3

# Could also define the function like this
# def lorenz_equs(t, y):
#     dy = [
#         SIGMA * (y[1] - y[0]),
#         -y[0] * y[2] + R * y[0] - y[1],
#         y[0] * y[1] - B * y[2],
#     ]
#     return dy

lorenz = lambda t,y: [SIGMA * (y[1] - y[0]), -y[0] * y[2] + R * y[0] - y[1], y[0] * y[1] - B * y[2]]

fig = plt.figure()
ax = fig.gca(projection='3d')
sol = solve_ivp(lorenz, [0, 50], [1, 1, 1], rtol=1e-8, atol=1e-10)
ax.plot(sol.y[0], sol.y[1], sol.y[2])
plt.show()
