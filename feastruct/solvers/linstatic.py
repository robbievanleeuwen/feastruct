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

    def __init__(self, analysis, analysis_cases, solver_settings=None):
        """Inits the LinearStatic class.

        :param analysis: Analysis object to solve
        :type analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
        :param analysis_cases: List of analysis cases to solve
        :type analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
        :param solver_settings: Settings to use in the solver - if not supplied, the default
            settings are adopted
        :type solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
        """

        super().__init__(
            analysis=analysis, analysis_cases=analysis_cases, solver_settings=solver_settings)

    def solve(self):
        """Executes the linear static finite element solver and saves the relevant results."""

        if self.solver_settings.linear_static.time_info:
            print('\n-Starting the linear static solver...\n')

        # assign the global degree of freedom numbers
        if self.solver_settings.linear_static.time_info:
            str = '--Assigning the global degree of freedom numbers...'
            self.function_timer(str, self.assign_dofs)
        else:
            self.assign_dofs()

        # assemble the global stiffness matrix
        if self.solver_settings.linear_static.time_info:
            str = '--Assembling the global stiffness matrix...'
            (K, _) = self.function_timer(str, self.assemble_stiff_matrix)
        else:
            (K, _) = self.assemble_stiff_matrix()

        # loop through each analysis case
        for (i, analysis_case) in enumerate(self.analysis_cases):
            if self.solver_settings.linear_static.time_info:
                print('\n--Analysis case {0}:'.format(i))

            # assemble the external force vector
            if self.solver_settings.linear_static.time_info:
                str = '---Assembling the external force vector...'
                (f_ext, f_eq) = self.function_timer(str, self.assemble_fext, analysis_case)
            else:
                (f_ext, f_eq) = self.assemble_fext(analysis_case=analysis_case)

            # apply the equivalent nodal loads
            f_ext -= f_eq

            # apply the boundary conditions
            (K_mod, f_ext) = self.apply_bcs(K=K, f_ext=f_ext, analysis_case=analysis_case)

            # solve for the displacement vector
            if self.solver_settings.linear_static.solver_type == 'direct':
                solver_func = self.direct_solver
            elif self.solver_settings.linear_static.solver_type == 'cgs':
                solver_func = self.cgs_solver

            if self.solver_settings.linear_static.time_info:
                str = '---Solving for the displacement vector using the {0} solver...'.format(
                    self.solver_settings.linear_static.solver_type)
                u = self.function_timer(str, solver_func, K_mod, f_ext)
            else:
                u = solver_func(K_mod, f_ext)

            # save the displacements to the DoF objects inside the Node objects
            self.save_displacements(u=u, analysis_case=analysis_case)

            # calculate the reaction forces
            if self.solver_settings.linear_static.time_info:
                str = '---Calculating reactions...'
                self.function_timer(str, self.calculate_reactions, K, u, f_eq, analysis_case)
            else:
                self.calculate_reactions(K=K, u=u, f_eq=f_eq, analysis_case=analysis_case)

            # calculate the element stresses
            if self.solver_settings.linear_static.time_info:
                str = '---Calculating element stresses...'
                self.function_timer(str, self.calculate_stresses, analysis_case)
            else:
                self.calculate_stresses(analysis_case=analysis_case)
