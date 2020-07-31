"""
GEL MODEL
FEniCS script for simulating gel around VICells
NW - Spring 2020

Issues
-Boundary conditions are prescribed on the midpoint of faces (triangles) as opposed to the
 vertices (nodes)
-Poorly handled (hardcoded) i/o for simulation data - file names, resolutions, etc.
"""

import meshio
from scipy.spatial import distance_matrix
import numpy as np
import pandas as pd
import sys
import time
import json
import progressbar
from time import sleep
from dolfin import *


## Define objects and functions
class bc_nw(UserExpression):

    def __init__(self, mesh, cell2trans_dict, **kwargs):
        self.mesh = mesh
        self._cell2trans_dict = cell2trans_dict
        self.cell_record = []
        self.x_record = []
        self.error_log = []
        super().__init__(**kwargs)

    def value_shape(self):
        return (3,)

    def eval_cell(self, value, x, cell):
        try:
            value[0], value[1], value[2] = self._cell2trans_dict[cell.index]
        except KeyError:
            value[0], value[1], value[2] = (0, 0, 0)
            self.error_log.append(cell)
        self.cell_record.append(cell)
        self.x_record.append(x)


class Surface(SubDomain):

    def init_record(self):
        self.x_record = []

    def inside(self, x, on_boundary):
        self.x_record.append(x)
        return on_boundary


def create_surf_midpoints(surf_mesh):
    """
    Was never able to prescribed boundary condition deformation at nodes.
    This is a hacky work-around.
    Args:
        surf_mesh:

    Returns:
        the physical midpoints of the faces that make up surf_mesh
    """
    '''
    cell_dict = dict(surf_mesh.cells)
    midpoints = np.zeros(cell_dict['triangle'].shape)
    for idx, triangle in enumerate(cell_dict['triangle']):
        midpoints[idx] = surf_mesh.points[cell_dict['triangle'][idx]].mean(0)
    '''
    midpoints = []
    for cell in cytod_surf.cells:
        if cell.type == "triangle":
            if len(midpoints) == 0:
                midpoints = cell.data
            else:
                midpoints = np.vstack([midpoints, cell.data])
    points = surf_mesh.points
    midpoint_coords = np.zeros((midpoints.shape[0], 3))
    for idx, face in enumerate(midpoints):
        midpoint_coords[idx, :] = np.mean((points[face[0]],
                                    points[face[1]],
                                    points[face[2]]), axis=0)
        #print(f"the first point of the first fapoints {points[face[0]]}")
        #print(f"the second point of the first fapoints {points[face[1]]}")
        #print(f"the third point of the first fapoints {points[face[2]]}")
    return midpoint_coords


## Import data
date_str = "sphere2ellipse_round_3"
gel_str = "1"
#path = "data/" + date_str + "_G" + gel_str + "/"
decimation_factor = 10 ## Arbitrarily scale down the boundary conditions. Raw BC displacement is too high.
output_folder = date_str + "_G" + gel_str + ""
cytod_surf = meshio.read("larger_sphere_manual" + ".msh")

cytod_faces = []
for cell in cytod_surf.cells:
    if cell.type == "triangle":
        if len(cytod_faces) == 0:
            cytod_faces = cell.data
        else:
            cytod_faces = np.vstack([cytod_faces, cell.data])

cytod_vol = meshio.read("bounding_box_manual" + ".msh")
print(f"cytod volume points length: {len(cytod_vol.points)}")
cytod_vol_faces = []
for cell in cytod_vol.cells:
    if cell.type == "tetra":
        if len(cytod_vol_faces) == 0:
            cytod_vol_faces = cell.data
        else:
            cytod_vol_faces = np.vstack([cytod_vol_faces, cell.data])
print("cytod volume faces:")
print(cytod_vol_faces)
print(f"length of cytod volume faces is {len(cytod_vol_faces)}")
mesh = Mesh()
with XDMFFile("bounding_box_manual" + ".xdmf") as infile:
    infile.read(mesh)

x_bound = 149.5
y_bound = 149.5
z_bound = 120
x_res = 0.29286
y_res = 0.29286
z_res = 0.8



surf_mesh1_midpoints = create_surf_midpoints(cytod_surf)
vert_disp = pd.read_csv("displacements_only_manual.csv",
                        header=None).values
vert_loc = pd.read_csv("displacements_location.csv",
                        header=None).values
#print(f"length of vert_disp is {len(vert_disp)}")
vert_disp = vert_disp/decimation_factor
midpoint_disp = np.zeros((cytod_faces.shape[0], 3))
midpoint_loc = np.zeros((cytod_faces.shape[0], 3))
for idx, face in enumerate(cytod_faces):
    if face[0] < 1481:
        if face[1] < 1481:
            if face[2] < 1481:
                midpoint_disp[idx, :] = np.mean((vert_disp[face[0]],
                                        vert_disp[face[1]],
                                        vert_disp[face[2]]), axis=0)
                midpoint_loc[idx, :] = np.mean((vert_loc[face[0]],
                                        vert_loc[face[1]],
                                        vert_loc[face[2]]), axis=0)

#print(f"The 3 nodes of the first face in cytod_faces are {cytod_faces[0]}")
face_0 = cytod_faces[0,0]
face_1 = cytod_faces[0,1]
face_2 = cytod_faces[0,2]
x = (vert_loc[face_0,0] + vert_loc[face_1,0] + vert_loc[face_2, 0])/3
#print(f"the x location of the first trangle midpoint should be {x}")
#print(f"the location of the first triangle of midpoint_loc is {midpoint_loc[0,:]}")
#print(f"size of midpoint_loc is {len(midpoint_loc)}")

# Reorganize midpoint_disp according to the location of the cytod_faces
midpoint_disp_new = np.zeros(midpoint_disp.shape)
print(midpoint_disp_new)
for index, face in enumerate(surf_mesh1_midpoints):
    x, y, z = [face[0], face[1], face[2]]
    dist_mat = distance_matrix(np.array([[x, y, z]]), midpoint_loc)
    loc = np.argmin(dist_mat)
    midpoint_disp_new[loc] = [midpoint_disp[loc,0], midpoint_disp[loc,1], midpoint_disp[loc,2]]
print(f"this is the new midpoint displacement array {midpoint_disp_new}")
print(f"this is the old midpoint displacement array {midpoint_disp}")
#print(f"length of midpoint_disp is {len(midpoint_disp)}")
mesh_boundaries = np.vstack((cytod_vol.points.min(0), cytod_vol.points.max(0))).T


## Mark surface facets of mesh (2D) AND
## Revert domain marking of OUTER boundary of gel && create cell_idx -> transformation dict
mvc = MeshValueCollection("size_t", mesh, 2)
subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc)
print(f"number of faces in Mesh Function Sizet {subdomains.size()}")
subdomains.set_all(0)
domains = MeshFunction("size_t", mesh, 2)
print(f"number of faces in Mesh Function {domains.size()}")
surf = Surface()
surf.init_record()
surf.mark(domains, 1)
blah = domains.array()
count = 0
for i in blah:
    if i == 1:
        count = count + 1
print(f"domain has {count} faces on its boundary")

change_log = []
cell_idx_list = np.zeros(midpoint_disp_new.shape[0])
cells_on_surface = []
print(f"mesh boundaries: {mesh_boundaries}")
meh = 0
for index, face in enumerate(faces(mesh)):
    print(index)
    x, y, z = face.midpoint().array()

    if domains.array()[index] == 1:
        if np.isclose(x, mesh_boundaries[0, 0], atol=1) or np.isclose(x, mesh_boundaries[0, 1], atol=1):
            domains.array()[index] = 0
            change_log.append(index)
        elif np.isclose(y, mesh_boundaries[1, 0], atol=1) or np.isclose(y, mesh_boundaries[1, 1], atol=1):
            domains.array()[index] = 0
            change_log.append(index)
        elif np.isclose(z, mesh_boundaries[2, 0], atol=1) or np.isclose(z, mesh_boundaries[2, 1], atol=1):
            domains.array()[index] = 0
            change_log.append(index)
        else:
            meh = meh + 1
            dist_mat = distance_matrix(np.array([[x, y, z]]), surf_mesh1_midpoints)
            cell_idx_list[np.argmin(dist_mat)] = face.entities(3)[0]
            blah = face.entities(3)[0]
            if blah < 324846:
                if len(cells_on_surface) == 0:
                    cells_on_surface = [cytod_vol_faces[0, 0], cytod_vol_faces[0, 1], cytod_vol_faces[0, 2], cytod_vol_faces[0, 3]]
                else:
                    cells_on_surface = np.vstack([cells_on_surface, [cytod_vol_faces[blah, 0], cytod_vol_faces[blah, 1], cytod_vol_faces[blah, 2], cytod_vol_faces[blah, 3]]])
    '''
    close = np.isclose([x,y,z], surf_mesh1_midpoints, atol=1)
    for i in range(0, len(close)):
        if all(close[i,:] == [True, True, True]):
            #print(f"I am here for the {meh} time")
            meh = meh + 1
            domains.array()[index] = 1
            #dist_mat = distance_matrix(np.array([[x, y, z]]), surf_mesh1_midpoints)
            #cell_idx_list[np.argmin(dist_mat)] = face.entities(3)[0]
            cell_idx_list[i] = face.entities(3)[0]
            blah = face.entities(3)[0]
            if blah < 324846:
                if len(cells_on_surface) == 0:
                    cells_on_surface = [cytod_vol_faces[0, 0], cytod_vol_faces[0, 1], cytod_vol_faces[0, 2], cytod_vol_faces[0, 3]]
                else:
                    cells_on_surface = np.vstack([cells_on_surface, [cytod_vol_faces[blah, 0], cytod_vol_faces[blah, 1], cytod_vol_faces[blah, 2], cytod_vol_faces[blah, 3]]])
            break
    '''
blah = domains.array()
count = 0
for i in blah:
    if i == 1:
        count = count + 1
print(f"domain has {count} faces on its boundary")

print(f"cells_on_surface {cells_on_surface}")

# Create inner surface mesh

cells = {"tetra": cells_on_surface}
points = cytod_vol.points
meshio.write_points_cells("inner_surface.vtk",
                            points,
                            cells)
'''
count = 0
blah = 0
for i in range(0, len(cell_idx_list)):
    temp = cell_idx_list[i]
    temp = int(temp)
    if temp < 324846:
        face = cytod_vol_faces[temp,:]
        cells_on_surface[count,0] = face[0] 
        cells_on_surface[count,1] = face[1] 
        cells_on_surface[count,2] = face[2]
        cells_on_surface[count,3] = face[3]
        count = count + 1
    if cells_on_surface[count, 0] == 0:
        blah = blah + 1
print(f"the number of zeros is {blah}")
print(f"number of faces in sphere .xdmf mesh is {len(cells_on_surface)}")
print(f"all the faces on the inner surface: {cells_on_surface}")
'''
#cell_idx_list = cell_idx_list[0:2962]
#print(np.size(cell_idx_list))
#print(np.size(midpoint_disp))
midpoint_test = np.zeros((len(cell_idx_list), 3))
midpoint_test[:,:] = 1
cell2trans_dict = dict(zip(cell_idx_list,
                           midpoint_test))
boundary_func = bc_nw(mesh, cell2trans_dict)



# Gel boundary conditions
V = VectorFunctionSpace(mesh, "Lagrange", 1)
sF = FunctionSpace(mesh, "Lagrange", 1)
du = TrialFunction(V)
v = TestFunction(V)
u = Function(V)
zero = Constant((0.0, 0.0, 0.0))
bcs = []
sbd = []
sbd.append(CompiledSubDomain("near(x[0], side)", side = 0))
sbd.append(CompiledSubDomain("near(x[1], side)", side = 0))
sbd.append(CompiledSubDomain("near(x[2], side)", side = 0))
sbd.append(CompiledSubDomain("near(x[0], side)", side = x_bound))
sbd.append(CompiledSubDomain("near(x[1], side)", side = y_bound))
sbd.append(CompiledSubDomain("near(x[2], side)", side = z_bound))
[bcs.append((DirichletBC(V, zero, sub))) for sub in sbd]
bcs.append(None)
bcs[-1] = DirichletBC(V, boundary_func, domains, 1)


## Setting up simulation

dx = Measure('dx', domain=mesh, subdomain_data=subdomains, metadata={'quadrature_degree': 2})


total_start = time.time()

d = len(u)
I = Identity(d)
F = I + grad(u)
B = Constant((0.0, 0.0, 0.0))
T = Constant((0.0, 0.0, 0.0))
C = F.T * F
Ic = tr(C)
J = det(F)

# Material properties
shr0, nu0 = 1.0, 0.45
mu = 3
lmbda = 1

# Math
psi = (mu / 2) * (Ic - 3) - mu * ln(J) + (lmbda / 2) * (ln(J)) ** 2
Pi = psi * dx - dot(B, u) * dx - dot(T, u) * ds
F = derivative(Pi, u, v)
J = derivative(F, u, du)

## Solving
solve(F == 0, u, bcs, J=J, solver_parameters={"newton_solver": {"linear_solver": "mumps"}})
print()
print("Total Time: ", time.time() - total_start)


## Exporting Data
hdf5_file = HDF5File(mesh.mpi_comm(),
                     "output/" + output_folder + "/function_dump.h5",
                     "w")
hdf5_file.write(u, "/function")
hdf5_file.close()

file = File("output/" + output_folder + "/solution.pvd")
file << u
