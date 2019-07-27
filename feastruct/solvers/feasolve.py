import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as linalg


class Solver:
    """Parent class for a finite element solver.

    Provides a number of methods suitable for most finite element analysis solver types.

    :cvar analysis: Analysis object to solve
    :vartype analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
    :cvar analysis_cases: List of analysis cases to solve
    :vartype analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
    :cvar solver_settings: Settings to use in the solver
    :vartype solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
    :cvar int ndof: Number of degrees of freedom in the analysis
    """

    def __init__(self, analysis, analysis_cases, solver_settings):
        """Inits the Solver class.

        :param analysis: Analysis object to solve
        :type analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
        :param analysis_cases: List of analysis cases to solve
        :type analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
        :param solver_settings: Settings to use in the solver
        :type solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
        """

        self.analysis = analysis
        self.analysis_cases = analysis_cases
        self.solver_settings = solver_settings

        # calculate total number of global dofs
        self.ndof = len(self.analysis.nodes) * len(self.analysis.dofs)

    def assign_dofs(self):
        """Method to assign global degrees of freedom to the nodes in the analysis object."""

        dof_count = 0  # initialise degree of freedom counter

        # loop through all the nodes in the analysis
        for node in self.analysis.nodes:
            # get the relevant degrees of freedom for the node
            node_dofs = node.get_dofs(self.analysis.dofs)

            # loop through each dof
            for dof in node_dofs:
                # assign a global degree of freedom number
                dof.global_dof_num = dof_count

                dof_count += 1

    def assemble_stiff_matrix(self, geometric=False, analysis_case=None):
        """Assembles the global stiffness using the sparse COO format.

        :param bool geometric: If set to True, also returns the global geometric stiffness matrix
        :param analysis_case: If geometric is set to True, the static analysis case from which to
            use to build the geometric stiffness matrix
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: The global stiffness matrix and, if desired, the global geometric stiffness
            matrix *(K, K_g)*
        :rtype: tuple(:class:`scipy.sparse.coo_matrix`, :class:`scipy.sparse.coo_matrix`)
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
            n = len(el.nodes) * len(self.analysis.dofs)

            # get element stiffness matrix
            k_el = el.get_stiffness_matrix()

            # get element degrees of freedom
            el_dofs = el.get_gdof_nums(dof_nums=self.analysis.dofs)

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
                k_el_g = el.get_geometric_stiff_matrix(analysis_case=analysis_case)
                k_g = k_el_g.flatten()
                data_g = np.hstack((data_g, k_g))

        K = sp.coo_matrix((data, (row, col)), shape=(self.ndof, self.ndof))

        if geometric:
            K_g = sp.coo_matrix((data_g, (row, col)), shape=(self.ndof, self.ndof))
        else:
            K_g = None

        return (K, K_g)

    def assemble_mass_matrix(self):
        """Assembles the global stiffness using the sparse COO format.

        :returns: The global mass matrix
        :rtype: :class:`scipy.sparse.coo_matrix`
        """

        # initialise lists
        row = []  # list containing row indices
        col = []  # list containing column indices
        data = []  # list containing mass matrix entries

        # loop through all the elements
        for el in self.analysis.elements:
            # determine number of dofs in the current element
            n = len(el.nodes) * len(self.analysis.dofs)

            # get element mass matrix
            m_el = el.get_mass_matrix()

            # get element degrees of freedom
            el_dofs = el.get_gdof_nums(dofs_nums=self.analysis.dofs)

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
        """Assembles the external force vector for analysis_case.

        :param analysis_case: Analysis case used to assemble the external force vector
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: The external force vector
        :rtype: :class:`numpy.ndarray`
        """

        f_ext = np.zeros(self.ndof)

        # apply nodal loads
        # loop through all the loads in the analysis case
        for load in analysis_case.load_case.items:
            load.apply_load(f_ext=f_ext)

        # add body forces
        for el in self.analysis.elements:
            # TODO: add body forces to fext
            pass

        return f_ext

    def apply_bcs(self, K, f_ext, analysis_case):
        """Applies the boundary conditions to the global stiffness matrix and external force
        vector for the analysis case.

        :param K: Global stiffness matrix
        :type K: :class:`scipy.sparse.coo_matrix`
        :param f_ext: External force vector
        :type f_ext: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: The global stiffness matrix and external force vector, modified by the boundary
            conditions *(K, f_ext)*
        :rtype: tuple(:class:`scipy.sparse.lil_matrix`, :class:`numpy.ndarray`)
        """

        # convert K to lil matrix
        K = sp.lil_matrix(K)

        # loop through all the supports in the analysis case
        for support in analysis_case.freedom_case.items:
            support.apply_support(K=K, f_ext=f_ext)

        # TODO: add spring stiffnesses

        return (K, f_ext)

    def remove_constrained_dofs(self, K, analysis_case):
        """Given an initial matrix *K*, returns a new matrix with the constrained degrees of
        freedom removed. The constrained degrees of freedom are those with supports defined in the
        freedom case used in the analysis case.

        :param K: Stiffness matrix
        :type K: :class:`scipy.sparse.coo_matrix`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: A new stiffness matrix with the constrined dofs removed
        :rtype: :class:`scipy.sparse.lil_matrix`
        """

        # convert K to lil matrix
        K_lil = sp.lil_matrix(K)

        # get dofs of constrained nodes
        constrained_dofs = []

        # loop through all the supports in the analysis case
        for support in analysis_case.freedom_case.items:
            constrained_dofs.append(support.node.dofs[support.dof].global_dof_num)

        # create list of free dofs
        free_dofs = [i for i in range(self.ndof) if i not in constrained_dofs]
        idx = np.ix_(free_dofs, free_dofs)

        return K_lil[idx]

    def direct_solver(self, K, f_ext):
        """Solves Ku=f_ext for *u* using the direct method.

        :param K: Global stiffness matrix
        :type K: :class:`scipy.sparse.lil_matrix`
        :param f_ext: External force vector
        :type f_ext: :class:`numpy.ndarray`

        :returns: Displacement vector *u*
        :rtype: :class:`numpy.ndarray`
        """

        # convert stiffness matrix to csc format
        K_csc = sp.csc_matrix(K)

        return linalg.spsolve(K_csc, f_ext)

    # def cgs_solver(self, K, f_ext):
    #     """Solves Ku=f_ext for *u* using the cgs method.
    #
    #     :param K: Global stiffness matrix
    #     :type K: :class:`scipy.sparse.lil_matrix`
    #     :param f_ext: External force vector
    #     :type f_ext: :class:`numpy.ndarray`
    #
    #     :returns: Displacement vector *u*
    #     :rtype: :class:`numpy.ndarray`
    #     """
    #
    #     # convert stiffness matrix to csc format
    #     K_csc = sp.csc_matrix(K)
    #
    #     # TODO: implement
    #     # if self.settings["precond"]:
    #     #     # perform ILU decomposition stiffness matrix
    #     #     pre_cond = linalg.LinearOperator(K.get_shape(),
    #     #                                      linalg.spilu(K_csc).solve)
    #     # else:
    #     #     pre_cond = None
    #     #
    #     # (u, exit) = linalg.cgs(K_csc, f_ext, tol=self.settings["tol"],
    #     #                        maxiter=self.settings["maxiter"], M=pre_cond)
    #     #
    #     # if (exit != 0):
    #     #     raise FEASolverError("CGS solver did not converge.")
    #
    #     # return u

    # def solve_eigenvalue(self, A, M):
    #     """
    #     """
    #
    #     A_csc = sp.csc_matrix(A)
    #     M_csc = sp.csc_matrix(M)
    #
    #     try:
    #         (w, v) = linalg.eigs(A=A_csc, k=self.settings["n"], M=M_csc,
    #                              sigma=self.settings["shift"],
    #                              maxiter=self.settings["maxiter"],
    #                              tol=self.settings["tol"])
    #
    #     except linalg.ArpackNoConvergence:
    #         raise FEASolverError("""Convergence not obtained for the
    #         eigenvalue solver""")
    #
    #     return (np.real(w), np.real(v))

    def save_displacements(self, u, analysis_case):
        """Saves the displacements *u* to the degree of freedom objects for the analysis case.

        :param u: Displacement vector
        :type u: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # loop through all the nodes in the analysis
        for node in self.analysis.nodes:
            # get the relevant degrees of freedom for the node
            node_dofs = node.get_dofs(self.analysis.dofs)

            # loop through each dof and save the relevant displacement
            for dof in node_dofs:
                dof.save_displacement(disp=u[dof.global_dof_num], analysis_case=analysis_case)

    # def save_eigenvectors(self, w, v, analysis_case, buckling=False,
    #                       frequency=False):
    #     """
    #     """
    #
    #     # add constrained dofs to eigenvector
    #     constrained_dofs = []
    #
    #     for support in analysis_case.freedom_case.items:
    #         constrained_dofs.append(support.node.dofs[support.dir-1])
    #
    #     for dof in sorted(constrained_dofs):
    #         v = np.insert(v, dof, 0, axis=0)
    #
    #     for node in self.analysis.nodes:
    #         for i in range(len(w)):
    #             w_i = w[i]
    #             v_i = v[node.dofs, i]
    #
    #             if buckling:
    #                 node.set_buckling_results(analysis_case.id, i+1, w_i, v_i)
    #             elif frequency:
    #                 node.set_frequency_results(analysis_case.id, i+1, w_i, v_i)

    def calculate_reactions(self, K, u, analysis_case):
        """Calculates the reactions using the stiffness matrix *K* and the displacement vector *u*
        and saves the reactions for the analysis case.

        :param K: Global stiffness matrix
        :type K: :class:`scipy.sparse.lil_matrix`
        :param u: Displacement vector
        :type u: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # calculate global force vector
        F = K.dot(u)

        # loop through constrained nodes and save reactions
        for support in analysis_case.freedom_case.items:
            f = F[support.node.dofs[support.dof].global_dof_num]

            support.save_reaction(f=f, analysis_case=analysis_case)

    def calculate_stresses(self, analysis_case):
        """Calculates and saves the element stresses for the analysis case.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # loop through each element in the analysis
        for el in self.analysis.elements:
            # get element stiffness matrix
            k_el = el.get_stiffness_matrix()

            # get nodal displacements
            u_el = el.get_nodal_displacements(analysis_case=analysis_case).reshape(-1)

            # calculate internal force vector
            f_int = np.matmul(k_el, u_el)

            el.save_fint(f=f_int, analysis_case=analysis_case)


class SolverSettings:
    """Class for settings to be used in a finite element analysis."""

    def __init__(self):
        """Inits the SolverSettings class."""

        pass
