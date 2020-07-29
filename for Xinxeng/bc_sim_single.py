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
    cell_dict = dict(surf_mesh.cells)
    midpoints = np.zeros(cell_dict['triangle'].shape)
    for idx, triangle in enumerate(cell_dict['triangle']):
        midpoints[idx] = surf_mesh.points[cell_dict['triangle'][idx]].mean(0)
    return midpoints


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
mesh = Mesh()
with XDMFFile("bounding_box_manual" + ".xdmf") as infile:
    infile.read(mesh)

mvc = MeshValueCollection("size_t", mesh, 2)

# Create derivative data
subdomains = cpp.mesh.MeshFunctionSizet(mesh, mvc)
subdomains.set_all(0)

x_bound = 149.5
y_bound = 149.5
z_bound = 120
x_res = 0.29286
y_res = 0.29286
z_res = 0.8

surf_mesh1_midpoints = create_surf_midpoints(cytod_surf)
vert_disp = pd.read_csv("displacements_only_manual.csv",
                        header=None).values

vert_disp = vert_disp/decimation_factor
midpoint_disp = np.zeros((cytod_faces.shape[0], 3))
for idx, face in enumerate(cytod_faces):
    midpoint_disp[idx, :] = np.mean((vert_disp[face[0]],
                                    vert_disp[face[1]],
                                    vert_disp[face[2]]), axis=0)
mesh_boundaries = np.vstack((cytod_vol.points.min(0), cytod_vol.points.max(0))).T
## Mark surface facets of mesh (2D)
domains = MeshFunction("size_t", mesh, 2)
surf = Surface()
surf.init_record()
surf.mark(domains, 1)


## Revert domain marking of OUTER boundary of gel && create cell_idx -> transformation dict
change_log = []
cell_idx_list = np.zeros(midpoint_disp.shape[0])
count = 0
for index, face in enumerate(faces(mesh)):
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
            count = count + 1
            dist_mat = distance_matrix(np.array([[x, y, z]]), surf_mesh1_midpoints)
            cell_idx_list[np.argmin(dist_mat)] = face.entities(3)[0]

## Setting up simulation

count = 0
for blah in cell_idx_list:
    if blah == 0:
        count = count + 1

dx = Measure('dx', domain=mesh, subdomain_data=subdomains, metadata={'quadrature_degree': 2})

V = VectorFunctionSpace(mesh, "Lagrange", 1)
sF = FunctionSpace(mesh, "Lagrange", 1)
du = TrialFunction(V)
v = TestFunction(V)
u = Function(V)

# Gel boundary conditions
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

total_start = time.time()
## Create boundary condition function
print(len(cell_idx_list))
print(len(midpoint_disp))
cell2trans_dict = dict(zip(cell_idx_list,
                           midpoint_disp))
boundary_func = bc_nw(mesh, cell2trans_dict)
bcs[-1] = DirichletBC(V, boundary_func, domains, 1)

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
