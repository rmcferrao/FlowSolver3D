from dolfin import *

def set_brinkman_problem(
                mesh,
                subdomains,
                boundaries,
                subodmainMarkers = {
                    'sub1': 1,
                    'sub2': 2
                },
                boundaryMarkers = {
                    'inlet': 1,
                    'outlet': 2,
                },
                mu=0.001003,  # Water Viscosity [Pa.s]
                k=9.869233e-15, #permeability of porous matrix
                pin=Constant(1.0),  # Input Pressure [Pa]
                pout=Constant(-1.0),  # Output Pressure  [Pa
                f=Constant((0.0, 0.0, 0.0))  # external force
                ):

        V = VectorElement('P', mesh.ufl_cell(), 2)
        Q = FiniteElement('P', mesh.ufl_cell(), 1)
        Element = MixedElement([V, Q])

        W = FunctionSpace(mesh, Element)
        (u, p) = TrialFunctions(W)
        (v, q) = TestFunctions(W)

        # Markers to integral infinitesimal terms
        ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
        dx = Measure('dx', domain=mesh, subdomain_data=subdomains)

        # normal vector
        n = FacetNormal(mesh)

        # boundary conditions
        bc1 = DirichletBC(W.sub(0).sub(1),Constant(0.0),boundaries,3) #y direction
        bc2 = DirichletBC(W.sub(0).sub(2),Constant(0.0),boundaries,4) #z direction
        bcs = [bc1, bc2]

        # Variational terms
        aD = lambda k, u, v: ((mu / k) * inner(u, v))
        aS = lambda u, v: 2 * mu * inner((grad(u)), grad(v))
        b = lambda u, q: (div(u) * q)
        lf = lambda f, v: inner(f, v)

        # complete variational form
        aB = aS(u, v) * dx(subodmainMarkers['sub1']) + aD(k, u, v) * dx(subodmainMarkers['sub2']) - b(v, p) * dx - b(u, q) * dx

        l = lf(f, v) * dx - pin * dot(v, n) * ds(boundaryMarkers['inlet']) - pout * dot(v, n) * ds(boundaryMarkers['outlet'])

        A = assemble(aB)
        L = assemble(l)
        for bc in bcs:
            bc.apply(A,L)

        return A, L, W

def solve_linear_system(A, L, W):
    w = Function(W)
    solve(A, w.vector(), L, 'mumps', 'none')

    return w
