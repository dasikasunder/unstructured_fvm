""" Script for finding and plotting the exact solution"""


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

fig = plt.figure()
ax = fig.gca(projection='3d')


# Setup problem parameters

L = 1.0     # Domain size
u_0 = 1.0   # Convection speed
alpha = 1.0 # Diffusion coeffieint
N = 100     # No. of grid points

x = np.linspace(0, L, N)
y = np.linspace(0, L, N)

X, Y = np.meshgrid(x, y) # Create a mesh using meshgrid

# Define few additional variables for simplicity

pi = np.pi
r1 = (u_0)/(2.0*alpha) + np.sqrt((pi)**2 + ((u_0)/(2.0*alpha))**2)
r2 = (u_0)/(2.0*alpha) - np.sqrt((pi)**2 + ((u_0)/(2.0*alpha))**2)

# Exact solution

PHI = (np.sin(pi*Y))*((r2*np.exp(r1*X + r2*L) - r1*np.exp(r1*L + r2*X))/(r2*np.exp(r2*L) - r1*np.exp(r1*L)))

# Plot the surface.
surf = ax.plot_surface(X, Y, PHI, cmap=cm.coolwarm, linewidth=0, antialiased=False)

# Customize the axes.

ax.set_title('Steady state solution')
ax.set_xlabel("$x$", fontsize=12)
ax.set_ylabel("$y$", fontsize=12)
ax.set_zlabel("$\phi$", fontsize=12)

# Add a color bar which maps values to colors.

fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
