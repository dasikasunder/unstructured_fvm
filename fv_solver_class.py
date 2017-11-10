import numpy as np
import numpy.matlib
import unstructured_grid_classes as ugc
import field_classes as fc

class Fv_solver:
    """
    Solver class for finite volume method
    """
    
    def __init__(self, gamma, phi):
        self.gamma = gamma
        self.phi = phi
        self.A = np.matlib.zeros((self.phi.n_c, self.phi.n_c)) # Assembly matrix
        self.b = np.matlib.zeros((self.phi.n_c, 1)) # RHS vector
    
    def assemble_matrix(self):
        
        mesh = self.phi._ScalarField__mesh
        self.phi.apply_boundary_conditions()
        self.phi.compute_values_at_interior_edge_centers()
        phi_grad_edges = self.phi.interpolate_gradients_to_edges()
        
        Cell_list = mesh.get_cell_list()
        
        # Iterate over all cells
        
        for i in range(0, self.phi.n_c):
            
            Cell = Cell_list[i]
            Cell_index = Cell.get_global_cell_number()
            E = [None]*3
            [E[0], E[1], E[2]] = Cell.get_cell_edges()
            
            for j in range(0,3):
                
                if (isinstance(E[j], ugc.Edge)):
                    
                    gDiff = E[j].gDiff()
                    T_f = E[j].vector_tf()
                    Edge_index = E[j].get_global_edge_number() - 1
                    
                    if(E[j].is_left_cell(Cell_index)):
                        pass
                    else:
                        T_f = T_f*(-1.0)
                        
                    self.A[i,i] += self.gamma*gDiff 
                    
                    [LC, RC] = E[j].get_straddling_cells()
                    
                    if (LC.get_global_cell_number() != Cell_index):
                        n_F = LC.get_global_cell_number()-1
                    else:
                        n_F = RC.get_global_cell_number()-1
                        
                    self.A[i, n_F] = -self.gamma*gDiff
                    
                    term = phi_grad_edges[Edge_index]*self.gamma

                    self.b[i] += term.dot(T_f) 
                    
                if (isinstance(E[j], ugc.boundaryEdge)):
                    
                    b_id = E[j].get_boundary_id()
                    
                    if (b_id == 2 or b_id == 21):
                        t_b = E[j].vector_tb()
                        Edge_index = E[j].get_global_edge_number() - 1
                        self.A[i,i] += self.gamma*E[j].gDiff()
                        term1 = self.gamma*E[j].gDiff()*self.phi.edge_center_value[Edge_index]
                        term2 = (phi_grad_edges[Edge_index]*self.gamma).dot(t_b)
                        FluxVb = term1 + term2
                        self.b[i] += FluxVb
                        
                    if (b_id == 3):
                        pass
                        
    def solve_using_direct_matrix_inversion(self):
        Fv_solver.assemble_matrix(self)
        phi_new = np.linalg.solve(self.A, self.b)
        
        for i in range(0, self.phi.n_c):
            self.phi.cell_centroid_value[i] = phi_new[i]