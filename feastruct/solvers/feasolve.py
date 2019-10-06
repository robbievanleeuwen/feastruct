import time
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
    :cvar int ndof: Number of active degrees of freedom in the analysis
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

        if solver_settings is None:
            self.solver_settings = SolverSettings()
        else:
            self.solver_settings = solver_settings

    def assign_dofs(self):
        """Method to assign global degrees of freedom to the nodes in the analysis object."""

        # assign node freedom allocations by looping through all elements
        for element in self.analysis.elements:
            element.apply_nfa()

        dof_count = 0  # initialise degree of freedom counter

        # loop through all the nodes in the analysis
        for node in self.analysis.nodes:
            # get all degree of freedoms in the node freedom signature
            node_dofs = node.get_dofs(freedom_signature=node.nfs)

            # loop through each dof
            for dof in node_dofs:
                # assign a global degree of freedom number
                dof.global_dof_num = dof_count

                dof_count += 1

        # assign the total number of global dofs
        self.ndof = dof_count

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
            n = el.get_ndof()

            # get element stiffness matrix
            k_el = el.get_stiffness_matrix()

            # get element degrees of freedom
            el_dofs = el.get_gdof_nums()

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
            n = el.get_ndof()

            # get element mass matrix
            m_el = el.get_mass_matrix()

            # get element degrees of freedom
            el_dofs = el.get_gdof_nums()

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

        :returns: The external force vector *f_ext* and the equivalent nodal loads vector *f_eq*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        f_ext = np.zeros(self.ndof)
        f_eq = np.zeros(self.ndof)

        # apply nodal loads
        # loop through all the nodal loads in the analysis case
        for nodal_load in analysis_case.load_case.items:
            nodal_load.apply_load(f_ext=f_ext)

        # apply element loads
        # loop through all the element loads in the analysis case
        for element_load in analysis_case.load_case.element_items:
            element_load.apply_load(f_eq=f_eq)

        # add body forces
        for el in self.analysis.elements:
            # TODO: add body forces to fext
            pass

        return (f_ext, f_eq)

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
            constrained_dofs.append(support.get_gdof())

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

    def cgs_solver(self, K, f_ext):
        """Solves Ku=f_ext for *u* using the cgs method.

        :param K: Global stiffness matrix
        :type K: :class:`scipy.sparse.lil_matrix`
        :param f_ext: External force vector
        :type f_ext: :class:`numpy.ndarray`

        :returns: Displacement vector *u*
        :rtype: :class:`numpy.ndarray`

        :raises Exception: If the cgs solver does not converge or if the preconditioner fails
        """

        # convert stiffness matrix to csc format
        K_csc = sp.csc_matrix(K)

        if self.solver_settings.linear_static.cgs_precond:
            # perform ILU decomposition stiffness matrix
            try:
                precond = linalg.LinearOperator(K.get_shape(), linalg.spilu(K_csc).solve)
            except RuntimeError:
                raise Exception('Preconditioner could not be constructed.')
        else:
            precond = None

        tol = self.solver_settings.linear_static.cgs_tol
        maxiter = self.solver_settings.linear_static.cgs_maxiter

        (u, exit) = linalg.cgs(K_csc, f_ext, tol=tol, maxiter=maxiter, M=precond)

        if (exit != 0):
            raise Exception('CGS solver did not converge.')

        return u

    def solve_eigenvalue(self, A, M, eigen_settings):
        """Finds eigenvalues and eigenvectors for A * x[i] = w[i] * M * x[i].

        :param A: A matrix
        :type A: :class:`scipy.sparse.coo_matrix`
        :param M: M matrix
        :type M: :class:`scipy.sparse.coo_matrix`
        :param eigen_settings: A settings object for use with the eigenvalue solver
        :type eigen_settings: :class:`feastruct.solvers.feasolve.LinearBucklingSettings` or
            :class:`feastruct.solvers.feasolve.NaturalFrequencySettings`

        :returns: Eigenvalues and corresponding eigenvectors *(w, v)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)

        :raises Exception: If the eigenvalue solver does not converge
        """

        A_csc = sp.csc_matrix(A)
        M_csc = sp.csc_matrix(M)

        n = eigen_settings.num_modes
        shift = eigen_settings.shift
        maxiter = eigen_settings.maxiter
        tol = eigen_settings.tol

        try:
            (w, v) = linalg.eigs(A=A_csc, k=n, M=M_csc, sigma=shift, maxiter=maxiter, tol=tol)

        except linalg.ArpackNoConvergence:
            raise Exception('Convergence not obtained for the eigenvalue solver')

        return (np.real(w), np.real(v))

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
            node_dofs = node.get_dofs(freedom_signature=node.nfs)

            # loop through each dof and save the relevant displacement
            for dof in node_dofs:
                dof.save_displacement(disp=u[dof.global_dof_num], analysis_case=analysis_case)

    def save_buckling_results(self, w, v, analysis_case):
        """Saves the buckling results *w* and *v* to the degree of freedom objects for the analysis
        case.

        :param w: Eigenvalues
        :type w: :class:`numpy.ndarray`
        :param v: Eigenvectors
        :type v: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # add constrained dofs to the eigenvectors
        constrained_dofs = []

        for support in analysis_case.freedom_case.items:
            constrained_dofs.append(support.get_gdof())

        # set eigenvector value to zero for constrained dof
        for dof in sorted(constrained_dofs):
            v = np.insert(v, dof, 0, axis=0)

        # save buckling results for each dof
        # loop through all the nodes
        for node in self.analysis.nodes:
            # loop through all the relevant dofs
            for dof in node.get_dofs(freedom_signature=node.nfs):
                # get global degree of freedom number
                gdof = dof.global_dof_num

                # list of buckling modes
                buckling_modes = list(range(len(w)))

                # get eigenvectors for the current degree of freedom
                v_dof = v[gdof, :]

                # save buckling results
                dof.save_buckling_modes(
                    buckling_modes=buckling_modes, w=w, v=v_dof, analysis_case=analysis_case)

    def save_frequency_results(self, w, v, analysis_case):
        """Saves the frequency results *w* and *v* to the degree of freedom objects for the
        analysis case.

        :param w: Eigenvalues
        :type w: :class:`numpy.ndarray`
        :param v: Eigenvectors
        :type v: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # add constrained dofs to the eigenvectors
        constrained_dofs = []

        for support in analysis_case.freedom_case.items:
            constrained_dofs.append(support.get_gdof())

        # set eigenvector value to zero for constrained dof
        for dof in sorted(constrained_dofs):
            v = np.insert(v, dof, 0, axis=0)

        # save frequency results for each dof
        # loop through all the nodes
        for node in self.analysis.nodes:
            # loop through all the relevant dofs
            for dof in node.get_dofs(freedom_signature=node.nfs):
                # get global degree of freedom number
                gdof = dof.global_dof_num

                # list of buckling modes
                frequency_modes = list(range(len(w)))

                # get eigenvectors for the current degree of freedom
                v_dof = v[gdof, :]

                # save buckling results
                dof.save_frequency_modes(
                    frequency_modes=frequency_modes, w=w, v=v_dof, analysis_case=analysis_case)

    def calculate_reactions(self, K, u, f_eq, analysis_case):
        """Calculates the reactions using the stiffness matrix *K* and the displacement vector *u*
        and saves the reactions for the analysis case.

        :param K: Global stiffness matrix
        :type K: :class:`scipy.sparse.lil_matrix`
        :param u: Displacement vector
        :type u: :class:`numpy.ndarray`
        :param f_eq: Global equivalent nodal loads vector of size *N*
        :type f_eq: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # calculate global force vector
        F = K.dot(u) + f_eq

        # loop through constrained nodes and save reactions
        for support in analysis_case.freedom_case.items:
            f = F[support.get_gdof()]

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

            # find element loads
            for element_load in analysis_case.load_case.element_items:
                # if the current element has an applied element load
                if element_load.element is el:
                    # add nodal equivalent loads to f_int
                    f_int += element_load.nodal_equivalent_loads()

            el.save_fint(f=f_int, analysis_case=analysis_case)

    def function_timer(self, text, function, *args):
        """Displays the message *text* and returns the time taken for a function, with arguments
        *args*, to execute. The value returned by the timed function is also returned.

        :param string text: Message to display
        :param function: Function to time and execute
        :type function: function
        :param args: Function arguments
        :returns: Value returned from the function
        """

        start_time = time.time()

        if text != '':
            print(text)

        result = function(*args)

        if text != '':
            print('----completed in {0:.6f} seconds---'.format(time.time() - start_time))

        return result


class SolverSettings:
    """Class for settings to be used in a finite element analysis.

    :cvar linear_static: Linear static settings
    :vartype linear_static: :class:`feastruct.solvers.feasolve.LinearStaticSettings`
    :cvar linear_buckling: Linear buckling settings
    :vartype linear_buckling: :class:`feastruct.solvers.feasolve.LinearBucklingSettings`
    :cvar natural_frequency: Natural frequency settings
    :vartype natural_frequency: :class:`feastruct.solvers.feasolve.NaturalFrequencySettings`
    """

    def __init__(self):
        """Inits the SolverSettings class."""

        self.linear_static = LinearStaticSettings()
        self.linear_buckling = LinearBucklingSettings()
        self.natural_frequency = NaturalFrequencySettings()


class LinearStaticSettings:
    """Class for settings to be used in a linear static analysis.

    :cvar bool time_info: If set to True, a detailed description of the computation and the time
        cost is printed to the terminal
    :cvar string solver_type: Solver type used to solve Ku = f - 'direct' or 'cgs'
    :cvar float cgs_tol: Tolerance to achieve for the cgs solver
    :cvar int cgs_maxiter: Maximum number of iterations for the cgs solver
    :cvar bool cgs_precond: If set to True, calculates a preconditioner for K using an ILU
        decomposition
    """

    def __init__(self, time_info=False, solver_type='direct', cgs_tol=1e-5, cgs_maxiter=None,
                 cgs_precond=True):
        """Inits the LinearStaticSettings class.

        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal
        :param string solver_type: Solver type used to solve Ku = f - 'direct' or 'cgs'
        :param float cgs_tol: Tolerance to achieve for the cgs solver
        :param int cgs_maxiter: Maximum number of iterations for the cgs solver
        :param bool cgs_precond: If set to True, calculates a preconditioner for K using an ILU
            decomposition
        """

        self.time_info = time_info
        self.solver_type = solver_type
        self.cgs_tol = cgs_tol
        self.cgs_maxiter = cgs_maxiter
        self.cgs_precond = cgs_precond


class LinearBucklingSettings:
    """Class for settings to be used in a linear buckling analysis."""

    def __init__(self, time_info=False, num_modes=4, shift=0, maxiter=None, tol=1e-9):
        """Inits the LinearBucklingSettings class.

        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal
        :param int num_modes: Number of modes to compute
        :param float shift: Finds eigenvalues near *shift* using the shift-invert mode
        :param int maxiter: Maximum number of iterations allowed
        :param float tol: Relative accuracy for eigenvalues (stopping criterion)
        """

        self.time_info = time_info
        self.num_modes = num_modes
        self.shift = shift
        self.maxiter = maxiter
        self.tol = tol


class NaturalFrequencySettings:
    """Class for settings to be used in a natural frequency analysis."""

    def __init__(self, time_info=False, num_modes=4, shift=0, maxiter=None, tol=1e-9):
        """Inits the NaturalFrequencySettings class.

        :param bool time_info: If set to True, a detailed description of the computation and the
            time cost is printed to the terminal
        :param int num_modes: Number of modes to compute
        :param float shift: Finds eigenvalues near *shift* using the shift-invert mode
        :param int maxiter: Maximum number of iterations allowed
        :param float tol: Relative accuracy for eigenvalues (stopping criterion)
        """

        self.time_info = time_info
        self.num_modes = num_modes
        self.shift = shift
        self.maxiter = maxiter
        self.tol = tol
