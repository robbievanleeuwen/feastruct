import sys
from solvers.feasolve import Solver
from fea.exceptions import FEAInputError, FEASolverError


class LinearStatic(Solver):
    """asdkjasd
    """

    def __init__(self, analysis, case_ids, solver='direct',
                 settings={"tol": 1e-5, "maxiter": None, "precond": True}):
        Solver.__init__(self, analysis, case_ids, solver, settings)

    def solve(self):
        try:
            self.assign_dofs()
            K = self.assemble_matrix()

            for case in self.case_ids:
                # get analysis case
                analysis_case = self.analysis.find_analysis_case(case)

                f_ext = self.assemble_fext(analysis_case)
                (K, f_ext) = self.apply_bcs(K, f_ext, analysis_case)

                if self.solver == 'direct':
                    u = self.direct_solver(K, f_ext)
                elif self.solver == 'cgs':
                    u = self.cgs_solver(K, f_ext)

                self.save_results(u, analysis_case)

        except FEAInputError as error:
            print("Error in linear static solver. {}".format(error))
            sys.exit(1)
        except FEASolverError as error:
            print("Error in linear static solver. {}".format(error))
            sys.exit(1)
