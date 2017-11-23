import numpy as np
import numpy.matlib
import unstructured_grid_classes as ugc
import field_classes as fc
from fvc_class import fvc

class advection_diffusion_equation_solver:
    """
    Solver class for solving Advection-Diffusion Equation
    """
    
    def __init__(self, gamma, u_0, dt, phi):
        self.gamma = gamma
        self.u_0 = u_0
        self.dt = dt
        self.phi = phi
        self.volume = np.zeros(self.phi.n_c)
        
        mesh = self.phi._ScalarField__mesh
        
        Cell_list = mesh.get_cell_list()
        for i in range(0, self.phi.n_c):
            self.volume[i] = Cell_list[i].area()
        
        self.A = np.matlib.zeros((self.phi.n_c, self.phi.n_c)) # Assembly matrix
        self.b = np.matlib.zeros((self.phi.n_c, 1)) # RHS vector

    
    def solve_using_explicit_euler_method(self):
        
        L = self.dt*(-fvc.div(self.u_0, self.phi) + fvc.laplacian(self.gamma, self.phi))/self.volume
        self.phi.cell_centroid_value += L
          
    def solve_using_rk3_tvd_method(self):
        phi_1 = self.phi
        phi_2 = self.phi
        L1 = self.dt*(-fvc.div(self.u_0, self.phi) + fvc.laplacian(self.gamma, self.phi))/self.volume
        phi_1.cell_centroid_value = self.phi.cell_centroid_value + L1
        L2 = self.dt*(-fvc.div(self.u_0, phi_1) + fvc.laplacian(self.gamma, phi_1))/self.volume
        phi_2.cell_centroid_value = (3./4.)*(self.phi.cell_centroid_value) + (1./4.)*(phi_1.cell_centroid_value) + (1./4.)*L2
        L3 = self.dt*(-fvc.div(self.u_0, phi_2) + fvc.laplacian(self.gamma, phi_2))/self.volume    
        self.phi.cell_centroid_value = (1./3.)*self.phi.cell_centroid_value + (2./3.)*(phi_2.cell_centroid_value) + (2./3.)*L3    
            
    def solve_using_rk4_method(self):
        phi_1 = self.phi
        phi_2 = self.phi
        phi_3 = self.phi

        L1 = self.dt*(-fvc.div(self.u_0, self.phi) + fvc.laplacian(self.gamma, self.phi))/self.volume

        phi_1.cell_centroid_value  = self.phi.cell_centroid_value + (1./4.)*L1

        L2 = self.dt*(-fvc.div(self.u_0, phi_1) + fvc.laplacian(self.gamma, phi_1))/self.volume

        phi_2.cell_centroid_value  = self.phi.cell_centroid_value + (1./3.)*L2
        
        L3 = self.dt*(-fvc.div(self.u_0, phi_2) + fvc.laplacian(self.gamma, phi_2))/self.volume
        
        phi_3.cell_centroid_value  = self.phi.cell_centroid_value + (1./2.)*L3
        
        L4 = self.dt*(-fvc.div(self.u_0, phi_3) + fvc.laplacian(self.gamma, phi_3))/self.volume
        
        self.phi.cell_centroid_value = self.phi.cell_centroid_value + L4
        
        
