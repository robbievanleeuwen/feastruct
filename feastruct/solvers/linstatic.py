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
            self.assign_dofs()  # assign the global degree of freedom numbers

            # assemble the global stiffness matrix
            (K, _) = self.assemble_stiff_matrix()

            # loop through each analysis case
            for case_id in self.case_ids:
                # get analysis case
                analysis_case = self.analysis.find_analysis_case(case_id)

                # assemble the external force vector
                f_ext = self.assemble_fext(analysis_case)

                # apply the boundary conditions
                (K_mod, f_ext) = self.apply_bcs(K, f_ext, analysis_case)

                # solve for the displacement vector
                if self.solver == 'direct':
                    u = self.direct_solver(K_mod, f_ext)
                elif self.solver == 'cgs':
                    u = self.cgs_solver(K_mod, f_ext)

                # save the displacements to the Node objects
                self.save_displacements(u, analysis_case)

                # calculate the reaction forces
                self.calculate_reactions(K, u, analysis_case)

                # calculate the element stresses
                self.calculate_stresses(analysis_case)

        except FEAInputError as error:
            print("Error in linear static solver. {}".format(error))
            sys.exit(1)
        except FEASolverError as error:
            print("Error in linear static solver. {}".format(error))
            sys.exit(1)
