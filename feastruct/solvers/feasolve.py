import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg
from feastruct.fea.exceptions import FEASolverError


class Solver:
    """adksa
    """

    def __init__(self, analysis, case_ids, solver, settings):
        """asdasdas
        """

        self.analysis = analysis
        self.case_ids = case_ids
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

    def assemble_stiff_matrix(self, geometric=False, case_id=None):
        """Assembles the global stiffness using the sparse COO format.
        """

        # initialise lists
        row = []  # list containing row indices
        col = []  # list containing column indices
        data = []  # list containing stiffness matrix entries

        if geometric:
            data_g = []

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

            # assemble geometric stiffness matrix
            if geometric:
                k_el_g = el.get_geometric_stiff_matrix(case_id)
                k_g = k_el_g.flatten()
                data_g = np.hstack((data_g, k_g))

        K = sp.coo_matrix((data, (row, col)), shape=(self.ndof, self.ndof))

        if geometric:
            K_g = sp.coo_matrix((data_g, (row, col)),
                                shape=(self.ndof, self.ndof))
        else:
            K_g = None

        return (K, K_g)

    def assemble_mass_matrix(self):
        """
        """

        # initialise lists
        row = []  # list containing row indices
        col = []  # list containing column indices
        data = []  # list containing mass matrix entries

        # loop through all the elements
        for el in self.analysis.elements:
            # determine number of dofs in the current element
            n = len(el.nodes) * self.analysis.dofs

            # get element mass matrix
            m_el = el.get_mass_matrix()

            # get element degrees of freedom
            el_dofs = el.get_dofs()

            # create row index vector
            r = np.repeat(el_dofs, n)

            # create column index vector
            c = np.tile(el_dofs, n)

            # flatten element mass matrix
            m = m_el.flatten()

            # add to global arrays
            row = np.hstack((row, r))
            col = np.hstack((col, c))
            data = np.hstack((data, m))

        return sp.coo_matrix((data, (row, col)), shape=(self.ndof, self.ndof))

    def assemble_fext(self, analysis_case):
        """asdsakd
        """

        f_ext = np.zeros(self.ndof)

        # apply loads
        for load in analysis_case.load_case.items:
            load.apply_load(f_ext)

        # add body forces
        for el in self.analysis.elements:
            # TODO: add body forces to fext
            pass

        return f_ext

    def apply_bcs(self, K, f_ext, analysis_case):
        """sdkljaskd

        expects K in coo_matrix format

        returns K in lil_matrix format
        """

        # convert K to lil matrix
        K_lil = sp.lil_matrix(K)

        for support in analysis_case.freedom_case.items:
            support.apply_support(K_lil, f_ext)

        # TODO: add spring stiffnesses

        return (K_lil, f_ext)

    def remove_constrained_dofs(self, K, analysis_case):
        """
        """

        # convert K to lil matrix
        K_lil = sp.lil_matrix(K)

        # get dofs of constrained nodes
        constrained_dofs = []

        for support in analysis_case.freedom_case.items:
            constrained_dofs.append(support.node.dofs[support.dir-1])

        # create list of free dofs
        free_dofs = [i for i in range(self.ndof) if i not in constrained_dofs]
        idx = np.ix_(free_dofs, free_dofs)

        return K_lil[idx]

    def direct_solver(self, K, f_ext):
        """asdkljaskdjsa

        expects K in lil_matrix format
        """

        # convert stiffness matrix to csc format
        K_csc = sp.csc_matrix(K)

        return linalg.spsolve(K_csc, f_ext)

    def cgs_solver(self, K, f_ext):
        """sadasdsa

        expects K in lil_matrix format
        """

        # convert stiffness matrix to csc format
        K_csc = sp.csc_matrix(K)

        if self.settings["precond"]:
            # perform ILU decomposition stiffness matrix
            pre_cond = linalg.LinearOperator(K.get_shape(),
                                             linalg.spilu(K_csc).solve)
        else:
            pre_cond = None

        (u, exit) = linalg.cgs(K_csc, f_ext, tol=self.settings["tol"],
                               maxiter=self.settings["maxiter"], M=pre_cond)

        if (exit != 0):
            raise FEASolverError("CGS solver did not converge.")

        return u

    def solve_eigenvalue(self, A, M):
        """
        """

        A_csc = sp.csc_matrix(A)
        M_csc = sp.csc_matrix(M)

        try:
            (w, v) = linalg.eigs(A=A_csc, k=self.settings["n"], M=M_csc,
                                 sigma=self.settings["shift"],
                                 maxiter=self.settings["maxiter"],
                                 tol=self.settings["tol"])

        except linalg.ArpackNoConvergence:
            raise FEASolverError("""Convergence not obtained for the
            eigenvalue solver""")

        return (np.real(w), np.real(v))

    def save_displacements(self, u, analysis_case):
        """ aslkdjlksad
        """

        for node in self.analysis.nodes:
            node.set_displacements(analysis_case.id, u[node.dofs])

    def save_eigenvectors(self, w, v, analysis_case, buckling=False,
                          frequency=False):
        """
        """

        # add constrained dofs to eigenvector
        constrained_dofs = []

        for support in analysis_case.freedom_case.items:
            constrained_dofs.append(support.node.dofs[support.dir-1])

        for dof in sorted(constrained_dofs):
            v = np.insert(v, dof, 0, axis=0)

        for node in self.analysis.nodes:
            for i in range(len(w)):
                w_i = w[i]
                v_i = v[node.dofs, i]

                if buckling:
                    node.set_buckling_results(analysis_case.id, i+1, w_i, v_i)
                elif frequency:
                    node.set_frequency_results(analysis_case.id, i+1, w_i, v_i)

    def calculate_reactions(self, K, u, analysis_case):
        """
        """

        # calculate global force vector
        F = K.dot(u)

        # loop through constrained nodes and save reactions
        for support in analysis_case.freedom_case.items:
            f = F[(support.node.dofs[support.dir-1])]
            support.set_reaction(analysis_case.id, f)

    def calculate_stresses(self, analysis_case):
        """
        """

        # loop through each element
        for el in self.analysis.elements:
            # get element stiffness matrix
            k_el = el.get_stiff_matrix()

            # get nodal displacements
            u_el = el.get_nodal_displacements(analysis_case.id).reshape(-1)

            # calculate internal force vector
            f_int = np.matmul(k_el, u_el)

            el.set_fint(analysis_case.id, f_int)
