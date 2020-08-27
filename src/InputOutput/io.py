import os
from dolfin import HDF5File, Mesh, MeshFunction, MPI, File

comm = MPI.comm_world

def mesh_to_h5(filePath, mesh, subdomains, boundaries):
    hdf = HDF5File(mesh.mpi_comm(), filePath, "w")
    hdf.write(mesh, "/mesh")
    hdf.write(subdomains, "/subdomains")
    hdf.write(boundaries, '/boundaries')
    hdf.close()

def read_h5(h5Path):
    mesh = Mesh()
    hdf = HDF5File(mesh.mpi_comm(), h5Path, 'r')
    hdf.read(mesh, '/mesh', False)
    subdomains = MeshFunction('size_t', mesh, 3)
    boundaries = MeshFunction('size_t', mesh, 2)

    hdf.read(subdomains, '/subdomains')
    hdf.read(boundaries, '/boundaries')

    return mesh, subdomains, boundaries

def save_results_hdf5(filePath, subdomains, boundaries, U):
    hdf = HDF5File(comm, filePath, "w")
    hdf.write(subdomains, "/subdomains")
    hdf.write(boundaries, "/boundaries")
    hdf.write(U, '/solution_vector')
    hdf.close()

def save_meshfiles_wrapper(PROJECT_ROOT, subFolder, mesh, subdomains, boundaries):
    filePath = os.path.join(PROJECT_ROOT, 'src', 'MeshFiles', 'mesh.h5')
    mesh_to_h5(filePath, mesh, subdomains, boundaries)

    meshVisuPath = os.path.join(PROJECT_ROOT, 'src', 'MeshFiles', subFolder, 'Visualization', 'mesh.pvd')
    File(meshVisuPath) << mesh

    sudomainsVisuPath = os.path.join(PROJECT_ROOT, 'src', 'MeshFiles', subFolder, 'Visualization', 'subdomains.pvd')
    File(sudomainsVisuPath) << subdomains

    boundariesVisuPath = os.path.join(PROJECT_ROOT, 'src', 'MeshFiles', subFolder, 'Visualization', 'boundaries.pvd')
    File(boundariesVisuPath) << boundaries

def save_results_wrapper(PROJECT_ROOT, subFolder, u, p, w, subdomains, boundaries):
    velocityVisuPath = os.path.join(PROJECT_ROOT, 'src', 'Results', subFolder, 'Visualization', 'velocity.pvd')
    File(velocityVisuPath) << u

    pressureVisuPath = os.path.join(PROJECT_ROOT, 'src', 'Results', subFolder, 'Visualization', 'pressure.pvd')
    File(pressureVisuPath) << p

    resultsH5Path = os.path.join(PROJECT_ROOT, 'Results', subFolder, 'solution.h5')
    save_results_hdf5(resultsH5Path, subdomains, boundaries, w)


