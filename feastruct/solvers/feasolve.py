import numpy as np
import scipy.sparse as sp


class Solver:
    """adksa
    """

    def __init__(self, analysis, solver_type):
        """asdasdas
        """

        self.analysis = analysis
        self.solver_type = solver_type

        # calculate total number of global dofs
        self.ndof = len(self.analysis.nodes) * self.analysis.dofs

    def assign_dofs(self):
        """aslkdjaks
        """

        dofs = np.array(range(self.analysis.dofs))
        dof_count = 0  # initialise degree of freedom counter

        for node in self.analysis.nodes:
            node.dofs = dofs + dof_count
            dof_count += self.analysis.dofs

    def assemble_matrix(self):
        """Assembles the global stiffness using the sparse COO format and is
        returned in the sparse CSC format.
        """

        # initialise lists
        row = []  # list containing row indices
        col = []  # list containing column indices
        data = []  # list containing stiffness matrix entries

        # loop through all the elements
        for el in self.analysis.elements:
            # determine number of dofs in the current element
            n = len(el.nodes) * self.analysis.dofs

            # get element stiffness matrix
            k_el = el.get_stiff_matrix()

            # get element degrees of freedom
            el_dofs = el.get_dofs()

            # create row index vector
            r = np.repeat(el_dofs, n)

            # create column index vector
            c = np.tile(el_dofs, n)

            # flatten element stiffness matrix
            k = k_el.flatten()

            # add to global arrays
            row = np.hstack((row, r))
            col = np.hstack((col, c))
            data = np.hstack((data, k))

        return sp.coo_matrix((data, (row, col)), shape=(self.ndof, self.ndof))

    def assemble_fext(self):
        """asdsakd
        """
        pass

    def apply_bcs(self):
        """sdkljaskd
        """
        pass

    def direct_solver(self):
        """asdkljaskdjsa
        """
        pass

    def cgs_solver(self):
        """sadasdsa
        """
        pass
