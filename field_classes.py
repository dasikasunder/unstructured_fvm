import unstructured_grid_classes as ugc
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from matplotlib import cm


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
               
        
        
    # Plotting functions
    def plot_3d_surf(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        xi = np.linspace(min(self.x_c), max(self.x_c), 33)
        yi = np.linspace(min(self.y_c), max(self.y_c), 33)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata(self.x_c, self.y_c, self.cell_centroid_value, xi, yi, interp='linear')
       
        surf = ax.plot_surface(X, Y, Z, rstride=5, cstride=5, cmap=cm.jet, vmin=np.nanmin(Z), vmax=np.nanmax(Z))
        ax.set_title('Steady state solution')
        ax.set_xlabel("$x$", fontsize=12)
        ax.set_ylabel("$y$", fontsize=12)
        ax.set_zlabel("$\phi$", fontsize=12)
        ax.set_zlim3d(np.min(Z), np.max(Z))
        fig.colorbar(surf)
        plt.show()
        
        
    def plot_gnuplot_surface_values(self):
        data_array = np.column_stack([self.x_e, self.y_e, self.edge_center_value])
        np.savetxt('plot.txt', data_array)

    def plot_to_text_file(self):
        data_array = np.column_stack([self.x_c, self.y_c, self.cell_centroid_value])
        np.savetxt('plot.txt', data_array)
        print("Data written to plot.txt")
        
        
