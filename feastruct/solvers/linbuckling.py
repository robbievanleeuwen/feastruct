from feastruct.solvers.feasolve import Solver


class LinearBuckling(Solver):
    """Class for a linear buckling solver.

    :cvar analysis: Analysis object to solve
    :vartype analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
    :cvar analysis_cases: List of analysis cases to solve
    :vartype analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
    :cvar solver_settings: Settings to use in the solver
    :vartype solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
    :cvar int ndof: Number of degrees of freedom in the analysis
    """

    def __init__(self, analysis, analysis_cases, solver_settings=None):
        """Inits the LinearBuckling class.

        :param analysis: Analysis object to solve
        :type analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
        :param analysis_cases: List of analysis cases to solve
        :type analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
        :param solver_settings: Settings to use in the solver - if not supplied, the default
            settings are adopted
        :type solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
        """

        super().__init__(analysis, analysis_cases, solver_settings)

        # TODO: check that analysis_case results available

    def solve(self):
        """Executes the linear buckling finite element solver and saves the relevant results."""

        if self.solver_settings.linear_buckling.time_info:
            print('\n-Starting the linear buckling solver...')

        # assign the global degree of freedom numbers
        if self.solver_settings.linear_static.time_info:
            str = '--Assigning the global degree of freedom numbers...'
            self.function_timer(str, self.assign_dofs)
        else:
            self.assign_dofs()

        # loop through each analysis case
        for (i, analysis_case) in enumerate(self.analysis_cases):
            if self.solver_settings.linear_buckling.time_info:
                print('\n--Analysis case {0}:'.format(i))

            # assemble the global stiffness and geometric stiffness matrices
            if self.solver_settings.linear_buckling.time_info:
                str = '---Assembling the global stiffness and geometric stiffness matrices...'
                (K, K_g) = self.function_timer(
                    str, self.assemble_stiff_matrix, True, analysis_case)
            else:
                (K, K_g) = self.assemble_stiff_matrix(geometric=True, analysis_case=analysis_case)

            # apply the boundary conditions
            K_mod = self.remove_constrained_dofs(K=K, analysis_case=analysis_case)
            K_mod_g = self.remove_constrained_dofs(K=K_g, analysis_case=analysis_case)

            # solve for the eigenvalues
            if self.solver_settings.linear_buckling.time_info:
                str = '---Solving for eigenvalues and eigenvectors ({0} modes)...'.format(
                    self.solver_settings.linear_buckling.num_modes)
                (w, v) = self.function_timer(
                    str, self.solve_eigenvalue, K_mod, -K_mod_g,
                    self.solver_settings.linear_buckling)
            else:
                (w, v) = self.solve_eigenvalue(
                    A=K_mod, M=-K_mod_g, eigen_settings=self.solver_settings.linear_buckling)

            self.save_buckling_results(w=w, v=v, analysis_case=analysis_case)
