#!/usr/bin/env python3

import numpy as np

def get_vertex_information_from_mesh():
    
    meshData = open("mesh_files/424_elements.geo", "r")

    no_of_vertices = int(meshData.readline())
    vertex_array = np.zeros(shape=(no_of_vertices,2))

    for i in range (0, no_of_vertices):
        vertex = meshData.readline()
        vertex_coordinates = vertex.split()
        vertex_array[i,0] = float(vertex_coordinates[0])
        vertex_array[i,1] = float(vertex_coordinates[1])

    meshData.close()
    return no_of_vertices, vertex_array

def get_cell_information_from_mesh():
    
    meshData = open("mesh_files/106_elements.geo", "r")
    
    no_of_vertices = int(meshData.readline())
    for i in range (0, no_of_vertices):
        meshData.readline()
    
    no_of_cells = int(meshData.readline())
    cell_array = np.zeros(shape=(no_of_cells,6))
    
    for i in range (0, no_of_cells):
        cell = meshData.readline()
        cell_data = cell.split()
        cell_array[i,0] = int(cell_data[1])
        cell_array[i,1] = int(cell_data[2])
        cell_array[i,2] = int(cell_data[3])
        
    for i in range (0, no_of_cells):
        cell = meshData.readline()
        cell_data = cell.split()
        cell_array[i,3] = int(cell_data[1])
        cell_array[i,4] = int(cell_data[2])
        cell_array[i,5] = int(cell_data[3])
        
    return no_of_cells, cell_array
