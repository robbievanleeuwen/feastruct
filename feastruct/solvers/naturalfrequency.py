import numpy as np
from feastruct.solvers.feasolve import Solver


class NaturalFrequency(Solver):
    """Class for a natural frequency solver.

    :cvar analysis: Analysis object to solve
    :vartype analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
    :cvar analysis_cases: List of analysis cases to solve
    :vartype analysis_cases: list[:class:`~feastruct.fea.cases.AnalysisCase`]
    :cvar solver_settings: Settings to use in the solver
    :vartype solver_settings: :class:`~feastruct.solvers.feasolve.SolverSettings`
    :cvar int ndof: Number of degrees of freedom in the analysis
    """

    def __init__(self, analysis, analysis_cases, solver_settings=None):
        """Inits the NaturalFrequency class.

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
        """Executes the natural frequency finite element solver and saves the relevant results."""

        if self.solver_settings.natural_frequency.time_info:
            print('\n-Starting the natural frequency solver...')

        # assign the global degree of freedom numbers
        self.assign_dofs()

        # loop through each analysis case
        for (i, analysis_case) in enumerate(self.analysis_cases):
            if self.solver_settings.natural_frequency.time_info:
                print('\n--Analysis case {0}:'.format(i))

            # assemble the global stiffness matrix
            if self.solver_settings.natural_frequency.time_info:
                str = '---Assembling the global stiffness matrix...'
                (K, _) = self.function_timer(str, self.assemble_stiff_matrix)
            else:
                (K, _) = self.assemble_stiff_matrix()

            # assemble the global mass matrix
            if self.solver_settings.natural_frequency.time_info:
                str = '---Assembling the global mass matrix...'
                M = self.function_timer(str, self.assemble_mass_matrix)
            else:
                M = self.assemble_mass_matrix()

            # apply the boundary conditions
            K_mod = self.remove_constrained_dofs(K=K, analysis_case=analysis_case)
            M_mod = self.remove_constrained_dofs(K=M, analysis_case=analysis_case)

            # solve for the eigenvalues
            if self.solver_settings.natural_frequency.time_info:
                str = '---Solving for eigenvalues and eigenvectors ({0} modes)...'.format(
                    self.solver_settings.natural_frequency.num_modes)
                (w, v) = self.function_timer(
                    str, self.solve_eigenvalue, K_mod, M_mod,
                    self.solver_settings.natural_frequency)
            else:
                (w, v) = self.solve_eigenvalue(
                    A=K_mod, M=M_mod, eigen_settings=self.solver_settings.natural_frequency)

            # compute natural frequencies in Hz
            w = np.sqrt(w) / 2 / np.pi

            self.save_frequency_results(w=w, v=v, analysis_case=analysis_case)
