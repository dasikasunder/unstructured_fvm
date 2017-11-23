from unstructured_grid_classes import Triangulation, Vector
from field_classes import ScalarField
from finite_volume_solvers import advection_diffusion_equation_solver
from fvc_class import fvc
"""
Name     : Dasika Sunder
S.R. No. : 14062
Dept.    : Mechanical
Degree   : PhD

About:
Driver script for solving Advection-Diffusion Equation on [0,1] X [0,1]
unstructured mesh using finite volume method.

Equation:
d(phi)/dt + u_0 (d(phi)/dx) = gamma*((d^2(phi)/dx^2) + (d^2(phi)/dy^2))

Boundary conditions:
Left boundary: phi = sin(pi*y)
Top and bottom walls: phi = 0
Right boundary: d(phi)/dx = 0
"""

# Input parameters
mesh_file_name = "106_elements.geo"
u_0 = Vector(0.0,0)
Gamma = 1.000
dt = 1e-4
finalTime = 10.0
no_of_time_steps = int(finalTime/dt)

# Set up the problem
mesh = Triangulation(mesh_file_name)
phi = ScalarField(mesh)
solver = advection_diffusion_equation_solver(Gamma, u_0, dt, phi)

time_counter = 0.0

for i in range(0, no_of_time_steps):
    solver.solve_using_explicit_euler_method()
    time_counter = i*dt
    
    if (i % 100 == 0):
        print("time = ", round(time_counter,2))

print("time = ", round(time_counter,2))
print("Done!")


# Plot solution
phi.plot_to_text_file()




