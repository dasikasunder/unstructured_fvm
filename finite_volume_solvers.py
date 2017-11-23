import numpy as np
import numpy.matlib
import unstructured_grid_classes as ugc
import field_classes as fc
from fvc_class import fvc

class advection_diffusion_equation_solver:
    """
    Solver class for Advection-Diffusion Equation
    """
    
    def __init__(self, gamma, u_0, dt, phi):
        self.gamma = gamma
        self.u_0 = u_0
        self.dt = dt
        self.phi = phi
        
        self.A = np.matlib.zeros((self.phi.n_c, self.phi.n_c)) # Assembly matrix
        self.b = np.matlib.zeros((self.phi.n_c, 1)) # RHS vector


            
    def solve_using_explicit_euler_method(self):
        
        mesh = self.phi._ScalarField__mesh 
        Cell_list = mesh.get_cell_list()
        laplace = fvc.laplacian(self.gamma, self.phi)
        divergence = fvc.div(self.u_0, self.phi)
        extra_term = np.zeros(self.phi.n_c)
        
        for i in range(0, self.phi.n_c):
            Cell = Cell_list[i]
            volume = Cell.area()
            extra_term[i] = self.phi.cell_centroid_value[i]*volume - self.dt*divergence[i] + self.dt*laplace[i]
            self.phi.cell_centroid_value[i] = extra_term[i]/volume
        
            
                     
                     
                
            
                
        
