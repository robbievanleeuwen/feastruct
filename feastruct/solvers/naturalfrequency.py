import sys
import numpy as np
from feastruct.solvers.feasolve import Solver
from feastruct.fea.exceptions import FEAInputError, FEASolverError


class NaturalFrequency(Solver):
    """asdkjasd
    """

    def __init__(self, analysis, case_id, solver=None,
                 settings={"n": 4, "shift": 0, "maxiter": None,
                           "tol": 1e-6}):
        Solver.__init__(self, analysis, case_id, solver, settings)

        # TODO: check that case_id results available
        # TODO: implement multiple cases

    def solve(self):
        try:
            self.assign_dofs()  # assign the global degree of freedom numbers

            # assemble the global stiffness matrix
            (K, _) = self.assemble_stiff_matrix()

            # assemble the global mass matrix
            M = self.assemble_mass_matrix()

            # get analysis case
            analysis_case = self.analysis.find_analysis_case(self.case_ids)

            # apply the boundary conditions
            K_mod = self.remove_constrained_dofs(K, analysis_case)
            M_mod = self.remove_constrained_dofs(M, analysis_case)

            # solve for the eigenvalues
            (w, v) = self.solve_eigenvalue(K_mod, M_mod)

            # compute natural frequencies in Hz
            w = np.sqrt(w) / 2 / np.pi

            self.save_eigenvectors(w, v, analysis_case, frequency=True)

        except FEAInputError as error:
            print("Error in linear buckling solver. {}".format(error))
            sys.exit(1)
        except FEASolverError as error:
            print("Error in linear buckling solver. {}".format(error))
            sys.exit(1)
