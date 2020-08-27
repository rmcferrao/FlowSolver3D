from dolfin import *
import numpy as np

class MarkSubdomain:
    def __init__(self, coquina_array, mesh, n_iterations = 1):
        self.coquina_array = coquina_array
        self.mesh = mesh
        self.voxels_len = [
            coquina_array.shape[0] - 1,
            coquina_array.shape[1] - 1,
            coquina_array.shape[2] - 1

        ]
        self.n_iteration = n_iterations
        self.delete_index = self.image_indexes(Cell(mesh, 0), 'index')

    def planes(self, p1, p2, p3, p4):
        p1_p2 = p2 - p1
        p1_p3 = p3 - p1
        p1_p4 = p4 - p1
        p2_p3 = p3 - p2
        p3_p4 = p4 - p3
        n1 = np.cross(p1_p2, p1_p3)
        n2 = np.cross(p1_p2, p1_p4)
        n3 = np.cross(p1_p3, p1_p4)
        n4 = np.cross(p3_p4, p2_p3)
        return [[n1, p1], [n2, p1], [n3, p1], [n4, p3]]

    def point_in_plane(self, ps, p1, p2, p3, p4):  # retorna se um ponto esta em um dos planos do tetra:
        # (Ponto_teste, 4 vertices do tetra): -> boolean
        index_boundary = []
        n_and_p = self.planes(p1, p2, p3, p4)
        for i in range(ps.shape[0]):
            p = ps[i]
            for n, p_in in n_and_p:
                plane_eq = lambda x, y, z: n[0] * (x - p_in[0]) + n[1] * (y - p_in[1]) + n[2] * (z - p_in[2])
                if plane_eq(p[0], p[1], p[2]) == 0:
                    index_boundary.append(i)
                    break
        return index_boundary

    def recursive_sum(self, res, coords, loc, len):
        for i in range(len - 1):
            res = np.vstack((res, (coords[loc] + coords[i + 1]) / 2))
        len = len - 1
        loc = loc + 1
        if len == 0:
            return res
        else:
            return self.recursive_sum(res, coords, loc, len)

    def image_indexes(self, cell, index_or_coords):
        coords = []
        coords.append([cell.midpoint().x(), cell.midpoint().y(), cell.midpoint().z()])
        coords.append(cell.get_vertex_coordinates()[:3])
        coords.append(cell.get_vertex_coordinates()[3:6])
        coords.append(cell.get_vertex_coordinates()[6:9])
        coords.append(cell.get_vertex_coordinates()[9:12])

        coords1 = np.array(coords)
        if index_or_coords == 'coords' and self.n_iteration == 1:
            return coords1
        elif index_or_coords == 'index':
            index_1 = self.point_in_plane(coords1, coords1[1], coords1[2], coords1[3], coords1[4])
            return index_1
            # return np.delete(coords1, index_1, axis=0)  # return middle point

        coords2 = self.recursive_sum(coords1, coords1, 0, coords1.shape[0])
        if index_or_coords == 'coords' and self.n_iteration == 2:
            return coords2
        elif index_or_coords == 'index':
            index_2 = self.point_in_plane(coords2, coords2[1], coords2[2], coords2[3], coords2[4])
            return index_2

        coords3 = self.recursive_sum(coords2, coords2, 0, coords2.shape[0])
        if index_or_coords == 'coords' and self.n_iteration == 3:
            return coords3
        elif index_or_coords == 'index':
            index_3 = self.point_in_plane(coords3, coords3[1], coords3[2], coords3[3], coords3[4])
            return index_3

        coords4 = self.recursive_sum(coords3, coords3, 0, coords3.shape[0])
        if index_or_coords == 'coords' and self.n_iteration == 4:
            return coords4
        elif index_or_coords == 'index':
            index_4 = self.point_in_plane(coords4, coords4[1], coords4[2], coords4[3], coords4[4])
            return index_4

    def sum_coords_in(self, coords, coquina_array):
        count = 0
        for coord in coords:
            point_image = int(coord[0] * self.voxels_len[0]), int(coord[1] * self.voxels_len[1]), int(coord[2] * self.voxels_len[2])
            if coquina_array[point_image] == 0:
                count += 1
        if (count / coords.shape[0]) > 0.5:
            return 1
        else:
            return 0

    def mark_inside(self, subdomains, marker = 2):
        for c in cells(self.mesh):
            index = c.index()
            coords_in = self.image_indexes(c, 'coords')
            coords_in = np.delete(coords_in, self.delete_index, axis=0)
            color_cell = self.sum_coords_in(coords_in, self.coquina_array)
            if color_cell == 1:
                subdomains.array()[index] = marker

        return subdomains

# Sub domain for inflow (right)
class Outflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > 1.0 - DOLFIN_EPS and on_boundary


# Sub domain for outflow (left)
class Inflow(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS and on_boundary


# Sub domain for outflow (left)
class YBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return (x[1] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS) and on_boundary


class ZBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return (x[2] > 1.0 - DOLFIN_EPS or x[2] < DOLFIN_EPS) and on_boundary