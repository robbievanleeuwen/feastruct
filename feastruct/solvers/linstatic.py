from solvers.feasolve import Solver


class LinearStatic(Solver):
    """asdkjasd
    """

    def __init__(self, analysis, solver_type='direct'):
        Solver.__init__(self, analysis)

    def solve(self):
        self.assign_dofs()
        K = self.assemble_matrix()
