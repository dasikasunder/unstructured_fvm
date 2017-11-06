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
        self.edge_center_value = np.ones(self.n_e)
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
    
    def compute_gradient_at_cell_centroids(self):
        """
        List[Vector]: Returns the gradients at cell centroids of the mesh
        """
        pass


    def plot(self):
        data_array = np.column_stack([self.x_c, self.y_c, self.cell_centroid_value])
        np.savetxt('plot.txt', data_array)