#!/usr/bin/env python3

from math import sqrt
import numpy as np
from stray_file import get_vertex_information_from_mesh

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
    

class Cell:
    def __init__(self, Vertex_A, Vertex_B, Vertex_C):
        self.A = Vertex_A
        self.B = Vertex_B
        self.C = Vertex_C
        
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
        rx = (xA + xB + xC)/3
        ry = (yA + yB + yC)/3
        r = Point(rx,ry)
        return r
         

class Triangulation:
    
    def get_vertex_information(self):
        vertex_array = []
        [no_of_vertices, vertex_coordinates] = get_vertex_information_from_mesh()
        for i in range(0, no_of_vertices):
            A = Vertex(vertex_coordinates[i,0], vertex_coordinates[i,1], i+1)
            vertex_array.append(A)   
        return vertex_array
 
 
mesh = Triangulation() 
v = mesh.get_vertex_information()
b = v[5].get_global_vertex_number()
print(b)