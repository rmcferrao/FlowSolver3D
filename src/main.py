from mesher.structuredMesh import MarkSubdomain, Inflow, Outflow, YBoundary, ZBoundary
from InputOutput.io import save_meshfiles_wrapper, save_results_wrapper
from Simulator.simulator import set_brinkman_problem, solve_linear_system
from skimage import io
import os
from dolfin import UnitCubeMesh, MeshFunction

PROJECT_ROOT = os.path.dirname(__file__)
PROJECT_ROOT = PROJECT_ROOT[:PROJECT_ROOT.rfind('/')]

if __name__ == '__main__':
    nx_elements, ny_elements, nz_elements = 17, 17, 17

    imgFile = 'rock4.tif'
    subFolder = imgFile[:imgFile.rfind('.')]

    imgPath = os.path.join(PROJECT_ROOT, 'image3D', imgFile)
    coquina_array = io.imread(imgPath)

    mesh = UnitCubeMesh(nx_elements, ny_elements, nz_elements)

    subdomains = MeshFunction('size_t', mesh, 3)
    subdomains.set_all(1)

    boundaries = MeshFunction('size_t', mesh, 2)
    boundaries.set_all(0)

    macro_pore = MarkSubdomain(coquina_array, mesh)
    macro_pore.mark_inside(subdomains)

    inflow, outflow = Inflow(), Outflow()
    ybounds, zbounds = YBoundary(), ZBoundary()

    inflow.mark(boundaries, 1)
    outflow.mark(boundaries, 2)
    ybounds.mark(boundaries, 3)
    zbounds.mark(boundaries, 4)

    save_meshfiles_wrapper(PROJECT_ROOT, subFolder, mesh, subdomains, boundaries)
    
    A, L, W = set_brinkman_problem(mesh, subdomains, boundaries)

    w = solve_linear_system(A, L, W)

    u,p = w.split()
    u.rename("u", "velocity")
    p.rename("p", "pressure")

    save_results_wrapper(PROJECT_ROOT, subFolder, u, p, w, subdomains, boundaries)