from math import sqrt
import numpy as np

# 2D Vector class

class Vector:
    """
    Vector class
    """
    # Constructor
    def __init__(self, x = 0, y = 0):
        """
        Constructor: Takes x and y component of the vector.
                     Default: x = 0.0, y = 0.0
        """
        self.x = x
        self.y = y

    # Get norm of a vector
    def norm(self):
        """
        Float: Finds the 2-norm (magnitutde) of the vector
        """
        return sqrt(self.x**2 + self.y**2)

    # Operators overloaded
    def __add__(self,other):
        """
        Vector: '+' operator
        """
        x = self.x + other.x
        y = self.y + other.y
        return Vector(x,y)

    def __sub__(self,other):
        """
        Vector: '-' operator
        """
        x = self.x - other.x
        y = self.y - other.y
        return Vector(x,y)

    def __mul__(self, a):
        """
        Vector: '*' operator (for multiplyng with a scalar)
        """
        x = a*self.x
        y = a*self.y
        return Vector(x,y)

    # Dot product with other vector
    def dot(self, other):
        """
        Float: Returns the dot product of two vectors
        """
        return ((self.x)*(other.x) + (self.y)*(other.y))

    def normalize(self):
        """
        Void: Converts the given vector into a unit vector
        """
        a = self.x
        b = self.y
        self.x = self.x/(sqrt(a**2 + b**2))
        self.y = self.y/(sqrt(a**2 + b**2))

# Point class

class Point:
    """
    Point class
    """
    def __init__(self, x_coordinate, y_coordinate):
        """
        Constructor: Takes x and y coordinates of the point
        """
        self.x = x_coordinate
        self.y = y_coordinate

    def get_coordinates(self):
        """
        Tuple[float]: Returns the x and y coordinates of the point
        """
        return self.x, self.y

    def distance_from_origin(self):
        """
        Float: Returns the distance of the point from origin
        """
        d = sqrt(self.x**2 + self.y**2)
        return d
    
    def distance_from_point(self, P):
        [xP, yP] = P.get_coordinates()
        d = sqrt((self.x - xP)**2 + (self.y - yP)**2)
        return d

# Vertex class inherited from point class - adds global vertex number

class Vertex(Point):
    """
    Vertex class (inhertited from point class)
    """
    def __init__(self, x_coordinate, y_coordinate, vertex_number):
        Point.__init__(self, x_coordinate, y_coordinate)
        self.global_vertex_number = vertex_number

    def get_global_vertex_number(self):
        """
        Int: Returns the global index of the vertex
        """
        return self.global_vertex_number

    def print_vertex_information(self):
        """
        Void: Prints important data about the point
        """
        print("x-coordinate:", self.x)
        print("y-coordinate:", self.y)
        print("Vertex number:", self.global_vertex_number)


# Edge class (for interior edges)

class Edge:
    """
    Edge class
    """

    # Constructor

    def __init__(self, Vertex_A, Vertex_B, edge_number, boundary_id, L_cell, R_cell = 0):
        """
        Construtor: Edge(Vertex_A, Vertex_B, edge_number, boundary_id, L_cell, R_cell)
            where,
            Vertex A: First point of the edge
            Vertex B: Second point of the edge
            edge_number: Global edge number
            boundary_id: Boundary ID of the edge
            L_cell: Left cell of the edge
            R_cell: Right cell of the edge
        """
        self.A = Vertex_A
        self.B = Vertex_B
        self.LC = L_cell
        self.RC = R_cell
        self.b_id = boundary_id
        self.global_edge_number = edge_number

    # Getter methods

    def get_edge_vertices(self):
        """
        Tuple[Vertex]: Returns the the vertices which make the edge
        """
        return self.A, self.B

    def get_straddling_cells(self):
        """
        Tuple[Cell]: Returns the left and right straddling cells respectively
        """
        return self.LC, self.RC

    def get_boundary_id(self):
        """
        int: Returns the boundary id
        """
        return self.b_id

    def get_global_edge_number(self):
        """
        int: Returns global edge number as give in the mesh file
        """
        return self.global_edge_number

    # Methods to get geometrical information of the edge

    def length(self):
        """
        float: Returns the length of the edge
        """
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        length = sqrt((xB - xA)**2 + (yB - yA)**2)
        return length

    def mid_point(self):
        """
        Point: Returns the mid-point of the edge
        """
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        xm = 0.5*(xA + xB)
        ym = 0.5*(yA + yB)
        M =  Point(xm,ym)
        return M

    def outward_unit_normal(self):
        """
        Vector: Returns outward unit normal of the edge - from left cell to right cell
        """
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        a = xA - xB
        b = yA - yB
        nx = -b/sqrt(a**2 + b**2)
        ny = a/sqrt(a**2 + b**2)
        n_v = Vector(nx, ny)
        return n_v

    def surface_area_vector(self):
        """
        Vector: Returns surface area vector of the edge - from left cell to right cell
        """
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        a = xA - xB
        b = yA - yB
        n_v = Vector(-b, a)
        return n_v

    def unit_tangent_vector(self):
        """
        Vector: Returns unit vector along the edge from vertex A to vertex B
        """
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        a = xA - xB
        b = yA - yB
        tx = a/sqrt(a**2 + b**2)
        ty = b/sqrt(a**2 + b**2)
        t_v = Vector(tx, ty)
        return t_v
    
    def unit_vector_joining_straddling_cell_centroids(self):
        """
        Vector: Returns a unit vector joining the left cell centroid and the 
        right cell centroid 
        """
        PL = self.LC.centroid()
        PR = self.RC.centroid()
        [xL, yL] = PL.get_coordinates()
        [xR, yR] = PR.get_coordinates()
        ex = xR - xL
        ey = yR - yL
        Ef = Vector(ex, ey)
        Ef.normalize()
        return Ef    


    def vector_joining_straddling_cell_centroids(self):
        """
        Vector: Returns a vector joining the left cell centroid and the right
        cell centroid - orthogonal/normal correction approach
        """
        PL = self.LC.centroid()
        PR = self.RC.centroid()
        [xL, yL] = PL.get_coordinates()
        [xR, yR] = PR.get_coordinates()
        ex = xR - xL
        ey = yR - yL
        Ef = Vector(ex, ey)
        Ef.normalize()
        l = Edge.length(self)
        Ef.x = l*Ef.x # Normal/Orthogonal correction
        Ef.y = l*Ef.y
        return Ef

    def vector_tf(self):
        """
        Vector: Returns Tf vector such that Sf = Ef + Tf
        """
        Ef = Edge.vector_joining_straddling_cell_centroids(self)
        Sf = Edge.surface_area_vector(self)
        Tf = Sf - Ef
        return Tf
    
    def distance_between_straddling_cell_centroids(self):
        """
        Float: Returns distance between the left cell centroid and right cell
        centroid
        """
        PL = self.LC.centroid()
        PR = self.RC.centroid()
        [xL, yL] = PL.get_coordinates()
        [xR, yR] = PR.get_coordinates()
        d_CF = sqrt((xR - xL)**2 + (yR - yL)**2)
        return d_CF
    
    def g_c(self):
        """
        Float: Returns interpolation factor g_c
        """
        PL = self.LC.centroid()
        PR = self.RC.centroid()
        EM = Edge.mid_point(self)
        
        [xL, yL] = PL.get_coordinates()
        [xR, yR] = PR.get_coordinates()
        [xM, yM] = EM.get_coordinates()
        
        Ff = sqrt((xR - xM)**2 + (yR - yM)**2)
        Cf = sqrt((xL - xM)**2 + (yL - yM)**2)
        
        gc = Ff/(Ff + Cf)
        
        return gc
    
    def gDiff(self):
        """
        Float: Returns geometric diffusion factor for the edge
        """
        Ef = Edge.length(self)
        d_CF = Edge.distance_between_straddling_cell_centroids(self)
        return Ef/d_CF
    
    def is_left_cell(self, cell_index):
        """
        Boolean: Returns True if given cell index is left cell index of the edge 
        """
        
        if(cell_index ==  self.LC.get_global_cell_number()):
            return True
        else:
            return False
        
    def print_edge_information(self):
        """
        Prints important information about the vertex
        """
        print("Edge number:", self.global_edge_number)
        print("Edge type: Interior edge")
        print("Edge vertices:", self.A.get_global_vertex_number(),",", self.B.get_global_vertex_number())
        print("   Vertex", self.A.get_global_vertex_number(), ":", self.A.get_coordinates())
        print("   Vertex", self.B.get_global_vertex_number(), ":", self.B.get_coordinates())
        print("Left cell:", self.LC.get_global_cell_number(), ", Right cell:", self.RC.get_global_cell_number())
        print("Edge Length:", Edge.length(self))
        M = Edge.mid_point(self)
        print("Mid point:", M.get_coordinates())
        s_f = Edge.surface_area_vector(self)
        print("Sf:",  s_f.x , s_f.y, )
        e_f = Edge.vector_joining_straddling_cell_centroids(self)
        print("Ef:", e_f.x, e_f.y)
        t_f = Edge.vector_tf(self)
        print("Tf:" , t_f.x, t_f.y)
        g_f = Edge.g_c(self)
        print("g_c =", g_f)


class boundaryEdge():

        """
        boundaryEdge class - similar to Edge class with few differences
        """
        # Constructor

        def __init__(self, Vertex_A, Vertex_B, edge_number, boundary_id, L_cell):
            """
            Construtor: Edge(Vertex_A, Vertex_B, edge_number, boundary_id, L_cell)
                where,
                Vertex A: First point of the edge
                Vertex B: Second point of the edge
                edge_number: Global edge number
                boundary_id: Boundary ID of the edge
                L_cell: Left cell of the edge
            """

            self.A = Vertex_A
            self.B = Vertex_B
            self.LC = L_cell
            self.b_id = boundary_id
            self.global_edge_number = edge_number

        # Getter methods

        def get_edge_vertices(self):
            """
            Tuple[Vertex]: Returns the the vertices which
            make the edge
            """
            return self.A, self.B

        def get_straddling_cells(self):
            """
            Cell: Returns the left straddling cells
            """
            return self.LC

        def get_boundary_id(self):
            """
            int: Returns the boundary id
            """
            return self.b_id

        def get_global_edge_number(self):
            """
            int: Returns global edge number as give in the mesh file
            """
            return self.global_edge_number

        # Methods to get geometrical information of the edge

        def length(self):
            """
            float: Returns the length of the edge
            """
            [xA,yA] = self.A.get_coordinates()
            [xB,yB] = self.B.get_coordinates()
            length = sqrt((xB - xA)**2 + (yB - yA)**2)
            return length

        def mid_point(self):
            """
            Point: Returns the mid-point of the edge
            """
            [xA,yA] = self.A.get_coordinates()
            [xB,yB] = self.B.get_coordinates()
            xm = 0.5*(xA + xB)
            ym = 0.5*(yA + yB)
            M =  Point(xm,ym)
            return M

        def outward_unit_normal(self):
            """Vector: Returns outward unit normal of the edge
            - from left cell to right cell
            """
            [xA,yA] = self.A.get_coordinates()
            [xB,yB] = self.B.get_coordinates()
            a = xA - xB
            b = yA - yB
            nx = b/sqrt(a**2 + b**2)
            ny = -a/sqrt(a**2 + b**2)
            n_v = Vector(nx, ny)
            return n_v

        def surface_area_vector(self):
            """
            Vector: Returns surface area vector of the edge - from left cell to
            to outsdie domain
            """
            [xA,yA] = self.A.get_coordinates()
            [xB,yB] = self.B.get_coordinates()
            a = xA - xB
            b = yA - yB
            n_v = Vector(-b, a)
            return n_v

        def unit_tangent(self):
            """
            Vector: Returns unit vector along the edge from vertex A
            to vertex B
            """
            [xA,yA] = self.A.get_coordinates()
            [xB,yB] = self.B.get_coordinates()
            a = xA - xB
            b = yA - yB
            tx = a/sqrt(a**2 + b**2)
            ty = b/sqrt(a**2 + b**2)
            t_v = Vector(tx, ty)
            return t_v

        def vector_joining_centroid_to_edge_mid_point(self):
       
            """
            Vector: Returns a vector joining the left cell centroid and edge 
            mid-point
            """
       
            PL = self.LC.centroid()
            EM = Edge.mid_point(self)
            [xL, yL] = PL.get_coordinates()
            [xR, yR] = EM.get_coordinates()
            ex = xR - xL
            ey = yR - yL
            Ef = Vector(ex, ey)
            Ef.normalize()
            l = boundaryEdge.length(self)
            Ef.x = l*Ef.x # Normal/Orthogonal correction
            Ef.y = l*Ef.y
            return Ef        
        
        def unit_vector_joining_centroid_to_edge_mid_point(self):
       
            """
            Vector: Returns a unit vector joining the left cell centroid and edge 
            mid-point
            """
       
            PL = self.LC.centroid()
            EM = Edge.mid_point(self)
            [xL, yL] = PL.get_coordinates()
            [xR, yR] = EM.get_coordinates()
            ex = xR - xL
            ey = yR - yL
            Ef = Vector(ex, ey)
            Ef.normalize()
            return Ef
        
        def vector_tb(self):
            """
            Vector: Returns Tf vector such that Sf = Ef + Tf
            """
            Ef = boundaryEdge.vector_joining_centroid_to_edge_mid_point(self)
            Sf = boundaryEdge.surface_area_vector(self)
            Tf = Sf - Ef
            return Tf
        
        def distance_between_straddling_cell_and_edge_mid_point(self):
            """
            Float: Returns distance between the left cell centroid and right cell
            centroid
            """
            PL = self.LC.centroid()
            EM = Edge.mid_point(self)
            [xL, yL] = PL.get_coordinates()
            [xR, yR] = EM.get_coordinates()
            d_CB = sqrt((xL - xR)**2 + (yL - yR)**2)
            return d_CB
        
        def gDiff(self):
            Eb = boundaryEdge.length(self)
            d_Cb = boundaryEdge.distance_between_straddling_cell_and_edge_mid_point(self)
            return Eb/d_Cb
            
        def print_edge_information(self):
            """
            Prints important information about the vertex
            """
            print("Edge number:", self.global_edge_number)
            print("Edge type: Boundary edge")
            print("Edge vertices:", self.A.get_global_vertex_number(),",", self.B.get_global_vertex_number())
            print("   Vertex", self.A.get_global_vertex_number(), ":", self.A.get_coordinates())
            print("   Vertex", self.B.get_global_vertex_number(), ":", self.B.get_coordinates())
            print("Left cell:", self.LC.get_global_cell_number())
            print("Edge Length:", boundaryEdge.length(self))
            M = boundaryEdge.mid_point(self)
            print("Mid point:", M.get_coordinates())
            s_f = boundaryEdge.surface_area_vector(self)
            print("Sf:",  s_f.x , s_f.y, )
            e_b = boundaryEdge.vector_joining_centroid_to_edge_mid_point(self)
            print("Ef:", e_b.x, e_b.y)
            t_b = boundaryEdge.vector_tb(self)
            print("Tf:" , t_b.x, t_b.y)
            


class Cell:
    def __init__(self, Vertex_A, Vertex_B, Vertex_C, cell_id, Edge_A = 0, Edge_B = 0, Edge_C = 0):
        self.A = Vertex_A
        self.B = Vertex_B
        self.C = Vertex_C
        self.EA = Edge_A
        self.EB = Edge_B
        self.EC = Edge_C
        self.c_id = cell_id

    # Getter methods

    def get_cell_vertices(self):
        return self.A, self.B, self.C

    def get_cell_edges(self):
        return self.EA, self.EB, self.EC

    def get_global_cell_number(self):
        return self.c_id

    # Setter method (for edges only)

    def set_cell_edges(self, Edge_A, Edge_B, Edge_C):
        self.EA = Edge_A
        self.EB = Edge_B
        self.EC = Edge_C

    # Methods to get geometrical information of the edge

    def area(self):
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        [xC,yC] = self.C.get_coordinates()
        area = 0.5*abs(xA*yB + xB*yC + xC*yA - xB*yA - xC*yB -xA*yC)
        return area

    def centroid(self):
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        [xC,yC] = self.C.get_coordinates()
        rx = (xA + xB + xC)/3.0
        ry = (yA + yB + yC)/3.0
        r = Point(rx,ry)
        return r
    
    def print_cell_information(self):
        """"
        Prints Critical Information about the cell
        """
        
        print("Cell number:", self.c_id)
        print("Cell Vertices:", self.A.get_global_vertex_number(),",", self.B.get_global_vertex_number(),",",self.C.get_global_vertex_number())
        print("   Vertex", self.A.get_global_vertex_number(), ":", self.A.get_coordinates())
        print("   Vertex", self.B.get_global_vertex_number(), ":", self.B.get_coordinates())
        print("   Vertex", self.C.get_global_vertex_number(), ":", self.C.get_coordinates())
        print("Cell edges:", self.EA.get_global_edge_number(),",", self.EB.get_global_edge_number(),",",self.EC.get_global_edge_number())
        Area = Cell.area(self)
        print("Cell area:", Area)
        M = Cell.centroid(self)
        print("Cell centroid:", M.get_coordinates())
        
        
class read_mesh_data:

    # Constructor

    def __init__(self, s):

        filepath  = "mesh_files/" + s
        meshData = open(filepath, "r")
        self.no_of_vertices = int(meshData.readline())
        self.vertex_array = np.zeros(shape=(self.no_of_vertices,2))

        for i in range (0, self.no_of_vertices):
            vertex = meshData.readline()
            vertex_coordinates = vertex.split()
            self.vertex_array[i,0] = float(vertex_coordinates[0])
            self.vertex_array[i,1] = float(vertex_coordinates[1])

        self.no_of_cells = int(meshData.readline())
        self.cell_array = np.zeros(shape=(self.no_of_cells,6))

        for i in range (0, self.no_of_cells):
            cell = meshData.readline()
            cell_data = cell.split()
            self.cell_array[i,0] = int(cell_data[1])
            self.cell_array[i,1] = int(cell_data[2])
            self.cell_array[i,2] = int(cell_data[3])

        for i in range (0, self.no_of_cells):
            cell = meshData.readline()
            cell_data = cell.split()
            self.cell_array[i,3] = int(cell_data[1])
            self.cell_array[i,4] = int(cell_data[2])
            self.cell_array[i,5] = int(cell_data[3])

        self.no_of_edges = int(meshData.readline())
        self.edge_array = np.zeros(shape=(self.no_of_edges,5))

        for i in range (0, self.no_of_edges):
            edge = meshData.readline()
            edge_data = edge.split()
            self.edge_array[i,0] = int(edge_data[0])
            self.edge_array[i,1] = int(edge_data[1])
            self.edge_array[i,2] = int(edge_data[2])
            self.edge_array[i,3] = int(edge_data[3])
            self.edge_array[i,4] = int(edge_data[4])

        meshData.close()

    # Getter methods

    def get_no_of_vertices_in_mesh(self):
        return self.no_of_vertices

    def get_no_of_cells_in_mesh(self):
        return self.no_of_cells

    def get_no_of_edges_in_mesh(self):
        return self.no_of_edges

    def get_vertex_information_from_mesh(self):
        return self.vertex_array

    def get_cell_information_from_mesh(self):
        return self.cell_array

    def get_edge_information_from_mesh(self):
        return self.edge_array
    
    


class Triangulation:

    # Class constructor

    def __init__(self, s):

        meshData = read_mesh_data(s)

        vertex_information = meshData.get_vertex_information_from_mesh()
        cell_information = meshData.get_cell_information_from_mesh()
        edge_information = meshData.get_edge_information_from_mesh()

        self.no_of_vertices = meshData.get_no_of_vertices_in_mesh()
        self.no_of_cells = meshData.get_no_of_cells_in_mesh()
        self.no_of_edges = meshData.get_no_of_edges_in_mesh()

        self.vertex_array = []

        for i in range(0, self.no_of_vertices):
            A = Vertex(vertex_information[i,0], vertex_information[i,1], i+1)
            self.vertex_array.append(A)

        self.cell_array = []

        for i in range(0, self.no_of_cells):
            A = self.vertex_array[int(cell_information[i,0] - 1)]
            B = self.vertex_array[int(cell_information[i,1] - 1)]
            C = self.vertex_array[int(cell_information[i,2] - 1)]
            cell_id = i+1
            cell_i = Cell(A, B, C, cell_id)
            self.cell_array.append(cell_i)

        self.edge_array = []
        self.number_of_boundary_edges = 0

        for i in range(0, self.no_of_edges):
            A = self.vertex_array[int(edge_information[i,0] - 1)]
            B = self.vertex_array[int(edge_information[i,1] - 1)]
            LC = self.cell_array[int(edge_information[i,2] - 1)]

            if (int(edge_information[i,3] - 1) != -1):
                RC = self.cell_array[int(edge_information[i,3] - 1)]
                b_id = int(edge_information[i,4])
                e_id = i+1
                E = Edge(A, B, e_id, b_id, LC, RC)
                self.edge_array.append(E)

            elif(int(edge_information[i,3] - 1) == -1):
                self.number_of_boundary_edges += 1
                b_id = int(edge_information[i,4])
                e_id = i+1
                EB = boundaryEdge(A, B, e_id, b_id, LC)
                self.edge_array.append(EB)


        for i in range(0, self.no_of_cells):
            EA = self.edge_array[int(cell_information[i,3] - 1)]
            EB = self.edge_array[int(cell_information[i,4] - 1)]
            EC = self.edge_array[int(cell_information[i,5] - 1)]
            self.cell_array[i].set_cell_edges(EA, EB, EC)

    # Getter methods

    def get_no_of_vertices(self):
        return self.no_of_vertices

    def get_vertex_list(self):
        return self.vertex_array

    def get_no_of_edges(self):
        return self.no_of_edges

    def get_edge_list(self):
        return self.edge_array

    def get_no_of_cells(self):
            return self.no_of_cells

    def get_cell_list(self):
            return self.cell_array

    # Methods to provide additional information about mesh

    def no_of_boundary_cells(self):
        n = 0
        for i in range(0, self.no_of_cells):
            Cell = self.cell_array[i]
            [EA, EB, EC] = Cell.get_cell_edges()
            b_idA = EA.get_boundary_id()
            b_idB = EB.get_boundary_id()
            b_idC = EC.get_boundary_id()
            if (b_idA + b_idB + b_idC != 0):
                n = n+1
        return n

    def no_of_boundary_edges(self):
        return self.number_of_boundary_edges

    def no_of_boundary_patches(self):
        counter = 1
        b_id_old = self.edge_array[0].get_boundary_id()
        for i in range(1, self.no_of_edges):
            Edge = self.edge_array[i]
            b_id = Edge.get_boundary_id()
            if (b_id_old != b_id):
                counter += 1
            b_id_old = b_id
        return counter

    def edge_iterator(self, a):

        """
        Tuple[int, int]: Returns the begining and ending index of edge with
        boundary ID 'a'. Note that indexing starts from 0.
        """

        begin_counter = 0
        end_counter = 0
        b_id_list = []

        for i in range(0, self.no_of_edges):
            Edge = self.edge_array[i]
            b_id = Edge.get_boundary_id()
            b_id_list.append(b_id)

        for i in range(0, self.no_of_edges):
            if (b_id_list[i] == a):
                begin_counter = i
                break

        counter = begin_counter

        while (b_id_list[counter] == a ):

            end_counter += 1
            counter += 1

            if (counter == (len(b_id_list))):
                break


        end_counter += begin_counter
        return begin_counter, end_counter-1
