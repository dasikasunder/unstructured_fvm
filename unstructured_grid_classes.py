from math import sqrt
import numpy as np

class Point:
    def __init__(self, x_coordinate, y_coordinate):
        self.x = x_coordinate
        self.y = y_coordinate

    def get_coordinates(self):
        return self.x, self.y

    def distance_from_origin(self):
        d = sqrt(self.x**2 + self.y**2)
        return d


class Vertex(Point):
    def __init__(self, x_coordinate, y_coordinate, vertex_number):
        Point.__init__(self, x_coordinate, y_coordinate)
        self.global_vertex_number = vertex_number

    def get_global_vertex_number(self):
        return self.global_vertex_number

# Edge class

class Edge:

    # Constructor

    def __init__(self, Vertex_A, Vertex_B, L_cell, R_cell, boundary_id, edge_number):
        self.A = Vertex_A
        self.B = Vertex_B
        self.LC = L_cell
        self.RC = R_cell
        self.b_id = boundary_id
        self.global_edge_number = edge_number

    # Getter methods

    def get_edge_vertices(self):
        return self.A, self.B

    def get_straddling_cells(self):
        return self.LC, self.RC

    def get_boundary_id(self):
        return self.b_id

    def get_global_edge_number(self):
        return self.global_edge_number

    # Methods to get geometrical information of the edge

    def length(self):
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        length = sqrt((xB - xA)**2 + (yB - yA)**2)
        return length

    def mid_point(self):
        [xA,yA] = self.A.get_coordinates()
        [xB,yB] = self.B.get_coordinates()
        xm = 0.5*(xA + xB)
        ym = 0.5*(yA + yB)
        M =  Point(xm,ym)
        return M


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

    # Methods to get geometrical information of the edge

    def area_of_cell(self):
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

class read_mesh_data:

    # Constructor

    def __init__(self):

        meshData = open("mesh_files/106_elements.geo", "r")
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

    def __init__(self):

        meshData = read_mesh_data()

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



    # Getter methods

    def get_no_of_vertices(self):
        return self.no_of_vertices

    def get_vertex_array(self):
        return self.vertex_array

    def get_no_of_cells(self):
            return self.no_of_cells

    def get_cell_array(self):
            return self.cell_array







mesh = Triangulation()
c = mesh.get_cell_array()
Cell2 = c[1]
[A, B, C] = Cell2.get_cell_vertices()

print(A.get_coordinates())
