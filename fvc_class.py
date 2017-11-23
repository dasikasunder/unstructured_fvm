import numpy as np
import unstructured_grid_classes as ugc

class fvc:
    
    def laplacian(Gamma, phi):
        
        laplace = np.zeros(phi.n_c)
        
        mesh = phi._ScalarField__mesh         
        
        Cell_list = mesh.get_cell_list()
        phi.apply_boundary_conditions()
        
        for i in range(0, phi.n_c):
            
            Cell = Cell_list[i]
            Cell_index = Cell.get_global_cell_number() - 1
            no_of_edges_in_cell = 3
            E = [None]*3
            [E[0], E[1], E[2]] = Cell.get_cell_edges()
            
            for j in range(0,no_of_edges_in_cell):
                
                if (isinstance(E[j], ugc.Edge)):
                    
                    d = E[j].distance_between_straddling_cell_centroids()
                    phi_N = 0.0
                    phi_P = 0.0
                    [LC, RC] = E[j].get_straddling_cells()
                    area = E[j].length()
                    left_cell_index = LC.get_global_cell_number() - 1
                    right_cell_index = RC.get_global_cell_number() - 1
                    
                    if (left_cell_index == Cell_index):
                        phi_P = phi.cell_centroid_value[left_cell_index]
                        phi_N = phi.cell_centroid_value[right_cell_index]
                    elif (right_cell_index == Cell_index):
                        phi_P = phi.cell_centroid_value[right_cell_index]
                        phi_N = phi.cell_centroid_value[left_cell_index]
                    
                    laplace[i] += Gamma*(area/d)*(phi_N - phi_P) 
                    
                elif (isinstance(E[j], ugc.boundaryEdge)):
                    
                    d = E[j].distance_between_straddling_cell_and_edge_mid_point()
                    edge_index = E[j].get_global_edge_number() -1
                    area = E[j].length()
                    
                    if (E[j].get_boundary_id() == 2 or E[j].get_boundary_id() == 21):
                        phi_b = phi.edge_center_value[edge_index]
                        phi_P = phi.cell_centroid_value[Cell_index]
                        laplace[i] += Gamma*(area/d)*(phi_b - phi_P) 
                        
                    elif (E[j].get_boundary_id() == 3):
                        laplace[i] += 0
                        
                        
        return laplace
    
    
    
    def div(U, phi):
        
        divergence = np.zeros(phi.n_c)
        
        mesh = phi._ScalarField__mesh         
        
        Cell_list = mesh.get_cell_list()
        phi.apply_boundary_conditions()
        
        for i in range(0, phi.n_c):
            
            Cell = Cell_list[i]
            Cell_index = Cell.get_global_cell_number() - 1
            no_of_edges_in_cell = 3
            E = [None]*3
            [E[0], E[1], E[2]] = Cell.get_cell_edges()
            
            for j in range(0,no_of_edges_in_cell):
                
                if (isinstance(E[j], ugc.Edge)):
                    
                    Sf = E[j].surface_area_vector()  
                    phi_N = 0.0
                    phi_P = 0.0
                    phi_f = 0.0
                    
                    [LC, RC] = E[j].get_straddling_cells()
                    left_cell_index = LC.get_global_cell_number() - 1
                    right_cell_index = RC.get_global_cell_number() - 1
                    
                    if (left_cell_index == Cell_index):
                        phi_P = phi.cell_centroid_value[left_cell_index]
                        phi_N = phi.cell_centroid_value[right_cell_index]
                    elif (right_cell_index == Cell_index):
                        phi_P = phi.cell_centroid_value[right_cell_index]
                        phi_N = phi.cell_centroid_value[left_cell_index]
                        Sf = Sf*(-1.0)
                    
                    F = U.dot(Sf)
                    
                    if (F >= 0):
                        phi_f = phi_P
                    elif (F < 0):
                        phi_f = phi_N
                    
                    divergence[i] += F*phi_f
                    
                elif (isinstance(E[j], ugc.boundaryEdge)):
                    
                    Sf = E[j].surface_area_vector()  
                    F = U.dot(Sf)
                    edge_index = E[j].get_global_edge_number() -1
                    
                    if (E[j].get_boundary_id() == 2 or E[j].get_boundary_id() == 21):
                        phi_b = phi.edge_center_value[edge_index]
                        phi_P = phi.cell_centroid_value[Cell_index]
                        divergence[i] += F*phi_b 
                        
                    elif (E[j].get_boundary_id() == 3):
                        divergence[i] += F*phi.cell_centroid_value[Cell_index]
        
        return divergence

