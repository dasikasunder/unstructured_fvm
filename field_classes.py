import unstructured_grid_classes as ugc
import numpy as np

mesh = ugc.Triangulation("6784_elements.geo")

class ScalarField:
    """
    ScalarField class
    """

    def __init__(self, mesh):
        """
        Constructor: Input - Object 'mesh' of Triangulation class
        """
        # Get information from the mesh
        self.n_v = mesh.get_no_of_vertices()
        self.n_e = mesh.get_no_of_edges()
        self.n_c = mesh.get_no_of_cells()

        # Private mesh object
        self.__mesh = mesh

        self.vertex_value = np.zeros(self.n_v)
        self.edge_center_value = np.zeros(self.n_e)
        self.cell_centroid_value = np.zeros(self.n_c)
        self.x_v = np.zeros(self.n_v)
        self.y_v = np.zeros(self.n_v)

        V = mesh.get_vertex_list()
        E = mesh.get_edge_list()
        C = mesh.get_cell_list()

        for i in range(0, self.n_v):
            self.x_v[i] = V[i].x
            self.y_v[i] = V[i].y

        self.x_e = np.zeros(self.n_e)
        self.y_e = np.zeros(self.n_e)

        for i in range(0, self.n_e):
            M = E[i].mid_point()
            self.x_e[i] = M.x
            self.y_e[i] = M.y

        self.x_c = np.zeros(self.n_c)
        self.y_c = np.zeros(self.n_c)

        for i in range(0, self.n_c):
            M = C[i].centroid()
            self.x_c[i] = M.x
            self.y_c[i] = M.y


    def apply_boundary_conditions(self):
        """
        Apply the appropriate boundary conditions for the boundary faces
        """
        E = self.__mesh.get_edge_list()

        # Top and bottom wall Dirichlet bcs (boundary_id = 21)
 
        e21_iterator = self.__mesh.edge_iterator(21)

        self.edge_center_value[e21_iterator[0]:e21_iterator[1]+1] = 0.0 
        
        # Left Dirichlet bc (boundary_id = 2)
       
        e2_iterator = self.__mesh.edge_iterator(2)

        b = np.sin(np.pi*self.y_e[e2_iterator[0]:e2_iterator[1]+1])

        self.edge_center_value[e2_iterator[0]:e2_iterator[1]+1] \
        = b
        
        # Right Neumann bc (Zero flux, boundary_id = 3)
        
        e3_iterator = self.__mesh.edge_iterator(3)
        
        for i in range(e3_iterator[0], e3_iterator[1]+1):
            LC = E[i].get_straddling_cells()
            n = LC.get_global_cell_number() - 1
            self.edge_center_value[i] = self.cell_centroid_value[n]
            
        pass
    
    def compute_values_at_interior_edge_centers(self):
        """
        Calculate value of scalar field at interior edge centers of the mesh
        """
        # Compute vertex values
        numerator = np.zeros(self.n_v)
        denominator = np.zeros(self.n_v)
        
        Cell_list = self.__mesh.get_cell_list()
        
        for i in range(0, self.n_c):
            Cell = Cell_list[i]
            C = Cell.centroid()
            no_of_vertices_in_cell = 3
            V = [None]*3
            [V[0], V[1], V[2]] = Cell.get_cell_vertices()
            
            for j in range(0, no_of_vertices_in_cell):
                n = V[j].get_global_vertex_number() - 1
                d_nC = V[j].distance_from_point(C)
                numerator[n] += self.cell_centroid_value[i]/d_nC
                denominator[n] += 1.0/d_nC
                
        for i in range(0, self.n_v):
            self.vertex_value[i] = numerator[i]/denominator[i]
        
        # Compute the values at edge centers
        
        Edge_list = self.__mesh.get_edge_list()
        
        for i in range(0, self.n_e):
            Edge = Edge_list[i]
            
            if (isinstance(Edge, ugc.Edge)):
                [V1, V2] = Edge.get_edge_vertices()
                nv1 = V1.get_global_vertex_number() - 1
                nv2 = V2.get_global_vertex_number() - 1
                phi_nv1 = self.vertex_value[nv1]
                phi_nv2 = self.vertex_value[nv2]
                
                self.edge_center_value[i] = 0.5*(phi_nv1 + phi_nv2)
                
                
            elif (isinstance(Edge, ugc.boundaryEdge)):
                pass
            
    def compute_gradient_at_cell_centroid(self):
        """
        List[Vector]: Compute the gradients at the cell centroids
        """
        Cell_list = self.__mesh.get_cell_list()
        phi_grad = []
        
        for i in range(0, self.n_c):
            
            Cell = Cell_list[i]
            no_of_edges_in_cell = 3
            edge_index = [None]*3
            Sf = [None]*3
            E = [None]*3
            [E[0], E[1], E[2]] = Cell.get_cell_edges()
            
            for j in range(0,no_of_edges_in_cell):
                
                edge_index[j] = E[j].get_global_edge_number() - 1
                
                if (isinstance(E[j], ugc.Edge)):
                    
                    [LC, RC] = E[j].get_straddling_cells() 
                    Sf[j] = E[j].surface_area_vector()
                    
                    if (LC.get_global_cell_number() == Cell.get_global_cell_number()):
                        pass
                    elif (RC.get_global_cell_number() == Cell.get_global_cell_number()):
                        Sf[j] = Sf[j]*(-1.0)
                        
                elif (isinstance(E[j], ugc.boundaryEdge)):
                    Sf[j] = E[j].surface_area_vector()
        
            grad = Sf[0]*self.edge_center_value[edge_index[0]] + \
                   Sf[1]*self.edge_center_value[edge_index[1]] + \
                   Sf[2]*self.edge_center_value[edge_index[2]]
             
            grad  = grad*(1.0/Cell.area())
            phi_grad.append(grad)
        
        return phi_grad
    
    def interpolate_gradients_to_edges(self):
        
        Edge_list = self.__mesh.get_edge_list()
        phi_grad_edges = [None]*self.n_e
        phi_grad_cell_centroids = ScalarField.compute_gradient_at_cell_centroid(self)
        
        # Interior edges
        
        e0_iterator = self.__mesh.edge_iterator(0)
        
        for i in range(e0_iterator[0], e0_iterator[1]+1):
            Edge = Edge_list[i]
            gC = Edge.g_c()
            [LC, RC] = Edge.get_straddling_cells()
            nL = LC.get_global_cell_number() - 1
            nR = RC.get_global_cell_number() - 1
            phi_grad_L = phi_grad_cell_centroids[nL]
            phi_grad_R = phi_grad_cell_centroids[nR]
            phi_grad_edges[i] = phi_grad_L*gC + phi_grad_R*(1.0-gC) 
            
            # Correction 
            
            phi_L = self.cell_centroid_value[nL]
            phi_R = self.cell_centroid_value[nR]
            d_CF = Edge.distance_between_straddling_cell_centroids()
            term1 = (phi_R - phi_L)/d_CF
            e_cf = Edge.unit_vector_joining_straddling_cell_centroids()
            term2 = e_cf.dot(phi_grad_edges[i])
            phi_grad_edges[i] += e_cf*(term1 - term2)
         
        # Inlet boundary
        
        e2_iterator = self.__mesh.edge_iterator(2)
        
        for i in range(e2_iterator[0], e2_iterator[1]+1):
            Edge = Edge_list[i]
            LC = Edge.get_straddling_cells()
            nL = LC.get_global_cell_number() - 1
            phi_grad_L = phi_grad_cell_centroids[nL]
            phi_grad_edges[i] = phi_grad_L
            
            # Correction
            
            phi_L = self.cell_centroid_value[nL]
            phi_b = self.edge_center_value[i]
            d_Cb = Edge.distance_between_straddling_cell_and_edge_mid_point()
            term1 = (phi_b - phi_L)/d_Cb
            e_cb = Edge.unit_vector_joining_centroid_to_edge_mid_point()
            term2 = e_cb.dot(phi_grad_edges[i])
            phi_grad_edges[i] += e_cb*(term1 - term2)
            
            
            
            
        # Top wall and bottom wall
        
        e21_iterator = self.__mesh.edge_iterator(21)
        
        for i in range(e21_iterator[0], e21_iterator[1]+1):
            Edge = Edge_list[i]
            LC = Edge.get_straddling_cells()
            nL = LC.get_global_cell_number() - 1
            phi_grad_L = phi_grad_cell_centroids[nL]
            phi_grad_edges[i] = phi_grad_L
            
            # Correction
            
            phi_L = self.cell_centroid_value[nL]
            phi_b = self.edge_center_value[i]
            d_Cb = Edge.distance_between_straddling_cell_and_edge_mid_point()
            term1 = (phi_b - phi_L)/d_Cb
            e_cb = Edge.unit_vector_joining_centroid_to_edge_mid_point()
            term2 = e_cb.dot(phi_grad_edges[i])
            phi_grad_edges[i] += e_cb*(term1 - term2)
            
        # Outlet boundary (No need for gradient at this boundary)
         
        e3_iterator = self.__mesh.edge_iterator(3)
        
        for i in range(e3_iterator[0], e3_iterator[1]+1):
            V = ugc.Vector(0.0, 0.0) # Dummy value
            phi_grad_edges[i] = V
            
        return phi_grad_edges
        
        
        
    def plot(self):
        data_array = np.column_stack([self.x_c, self.y_c, self.cell_centroid_value])
        np.savetxt('plot.txt', data_array)