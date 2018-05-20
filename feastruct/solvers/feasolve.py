import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve, cgs, LinearOperator, spilu


class Solver:
    """adksa
    """

    def __init__(self, analysis, solver, settings):
        """asdasdas
        """

        self.analysis = analysis
        self.solver = solver
        self.settings = settings

        # calculate total number of global dofs
        self.ndof = len(self.analysis.nodes) * self.analysis.dofs

    def assign_dofs(self):
        """aslkdjaks
        """

        dofs = np.array(range(self.analysis.dofs))
        dof_count = 0  # initialise degree of freedom counter

        for node in self.analysis.nodes:
            node.dofs = dofs + dof_count
            dof_count += self.analysis.dofs

    def assemble_matrix(self):
        """Assembles the global stiffness using the sparse COO format.
        """

        # initialise lists
        row = []  # list containing row indices
        col = []  # list containing column indices
        data = []  # list containing stiffness matrix entries

        # loop through all the elements
        for el in self.analysis.elements:
            # determine number of dofs in the current element
            n = len(el.nodes) * self.analysis.dofs

            # get element stiffness matrix
            k_el = el.get_stiff_matrix()

            # get element degrees of freedom
            el_dofs = el.get_dofs()

            # create row index vector
            r = np.repeat(el_dofs, n)

            # create column index vector
            c = np.tile(el_dofs, n)

            # flatten element stiffness matrix
            k = k_el.flatten()

            # add to global arrays
            row = np.hstack((row, r))
            col = np.hstack((col, c))
            data = np.hstack((data, k))

        return sp.coo_matrix((data, (row, col)), shape=(self.ndof, self.ndof))

    def assemble_fext(self):
        """asdsakd
        """

        f_ext = np.zeros(self.ndof)

        # add nodal loads
        for nodal_load in self.analysis.nodal_loads:
            dofs = nodal_load["node"].dofs  # get dofs for current node
            dir = nodal_load["dir"]  # get direction of load

            # add load to f_ext, selecting the correct dof from dofs
            f_ext[dofs[dir-1]] = nodal_load["val"]

        # add body forces
        for el in self.analysis.elements:
            # TODO: add body forces to fext
            pass

        return f_ext

    def apply_bcs(self, K, f_ext):
        """sdkljaskd

        expects K in coo_matrix format

        returns K in lil_matrix format
        """

        # convert K to lil matrix
        K_lil = sp.lil_matrix(K)

        for support in self.analysis.supports:
            dofs = support["node"].dofs  # get dofs for current node
            dir = support["dir"]  # get direction of load

            # modify stiffness matrix and f_ext
            K_lil[dofs[dir-1], :] = 0
            K_lil[dofs[dir-1], dofs[dir-1]] = 1
            f_ext[dofs[dir-1]] = support["val"]

        # TODO: add spring stiffnesses

        return (K_lil, f_ext)

    def direct_solver(self, K, f_ext):
        """asdkljaskdjsa

        expects K in lil_matrix format
        """

        # convert stiffness matrix to csc format
        K_csc = sp.csc_matrix(K)

        return spsolve(K_csc, f_ext)

    def cgs_solver(self, K, f_ext):
        """sadasdsa

        expects K in lil_matrix format
        """

        # convert stiffness matrix to csc format
        K_csc = sp.csc_matrix(K)

        if self.settings["precond"]:
            # perform ILU decomposition stiffness matrix
            pre_cond = LinearOperator(K.get_shape(), spilu(K_csc).solve)
        else:
            pre_cond = None

        (u, exit) = cgs(K_csc, f_ext, tol=self.settings["tol"],
                        maxiter=self.settings["maxiter"], M=pre_cond)

        if (exit != 0):
            raise RuntimeError("CGS solver did not converge.")

        return u

    def save_results(self, u):
        """ aslkdjlksad
        """

        for node in self.analysis.nodes:
            node.u = u[node.dofs]
