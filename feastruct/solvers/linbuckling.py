import sys
from feastruct.solvers.feasolve import Solver
from feastruct.fea.exceptions import FEAInputError, FEASolverError


class LinearBuckling(Solver):
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
            # assemble the global stiffness and geometric stiffnes matrices
            (K, K_g) = self.assemble_stiff_matrix(geometric=True,
                                                  case_id=self.case_ids)

            # get analysis case
            analysis_case = self.analysis.find_analysis_case(self.case_ids)

            # apply the boundary conditions
            K_mod = self.remove_constrained_dofs(K, analysis_case)
            K_mod_g = self.remove_constrained_dofs(K_g, analysis_case)

            # solve for the eigenvalues
            (w, v) = self.solve_eigenvalue(K_mod, -K_mod_g)

            self.save_eigenvectors(w, v, analysis_case, buckling=True)

        except FEAInputError as error:
            print("Error in linear buckling solver. {}".format(error))
            sys.exit(1)
        except FEASolverError as error:
            print("Error in linear buckling solver. {}".format(error))
            sys.exit(1)
