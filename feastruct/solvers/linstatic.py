import sys
from solvers.feasolve import Solver


class LinearStatic(Solver):
    """asdkjasd
    """

    def __init__(self, analysis, solver='direct',
                 settings={"tol": 1e-5, "maxiter": None, "precond": True}):
        Solver.__init__(self, analysis, solver, settings)

    def solve(self):
        try:
            self.assign_dofs()
            K = self.assemble_matrix()
            f_ext = self.assemble_fext()
            (K, f_ext) = self.apply_bcs(K, f_ext)

            if self.solver == 'direct':
                u = self.direct_solver(K, f_ext)
            elif self.solver == 'cgs':
                u = self.cgs_solver(K, f_ext)

        except RuntimeError as error:
            print("Error in linear static solver. {}".format(error))
            sys.exit(1)

        self.save_results(u)
        print(u)
