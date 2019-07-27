from feastruct.solvers.feasolve import Solver


class LinearStatic(Solver):
    """Class for a linear static solver.

    :cvar analysis: Analysis object to solve
    :vartype analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
    :cvar analysis_cases: List of analysis cases to solve
    :vartype analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
    :cvar solver_settings: Settings to use in the solver
    :vartype solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
    :cvar int ndof: Number of degrees of freedom in the analysis
    """

    def __init__(self, analysis, analysis_cases, solver_settings):
        """Inits the LinearStatic class.

        :param analysis: Analysis object to solve
        :type analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
        :param analysis_cases: List of analysis cases to solve
        :type analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
        :param solver_settings: Settings to use in the solver
        :type solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
        """

        super().__init__(
            analysis=analysis, analysis_cases=analysis_cases, solver_settings=solver_settings)

    def solve(self):
        """Executes the finite element solver and saves the relevant results."""

        # assign the global degree of freedom numbers
        self.assign_dofs()

        # assemble the global stiffness matrix
        (K, _) = self.assemble_stiff_matrix()

        # loop through each analysis case
        for analysis_case in self.analysis_cases:
            # assemble the external force vector
            f_ext = self.assemble_fext(analysis_case=analysis_case)

            # apply the boundary conditions
            (K_mod, f_ext) = self.apply_bcs(K=K, f_ext=f_ext, analysis_case=analysis_case)

            # solve for the displacement vector
            u = self.direct_solver(K=K_mod, f_ext=f_ext)
            # TODO: implement cgs solver
            # elif self.solver == 'cgs':
            #     u = self.cgs_solver(K_mod, f_ext)

            # save the displacements to the DoF objects inside the Node objects
            self.save_displacements(u=u, analysis_case=analysis_case)

            # calculate the reaction forces
            self.calculate_reactions(K=K, u=u, analysis_case=analysis_case)

            # calculate the element stresses
            self.calculate_stresses(analysis_case=analysis_case)
