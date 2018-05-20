import numpy as np
from fea.fea import fea, FiniteElement


class Frame2D(fea):
    """
    """

    def __init__(self, analysis_type='linear', nodes=[], elements=[]):
        """
        """

        # initialise parent fea class
        fea.__init__(self, analysis_type, nodes, elements)
        self.dofs = 3

    def add_element(self, id, node_ids, el_type='EB2', E=1, A=1, ixx=1,
                    G=0, A_s=0):
        """
        EB2: 2 Noded Euler Bernoulli Beam
        TB2: 2 Noded Timoshenko Beam
        """

        # TODO: check value types e.g. id and node_ids are ints

        if el_type == 'EB2':
            element = EulerBernoulliBeam2D(self, id, node_ids, E, A, ixx)
        elif el_type == 'TB2':
            element = TimoshenkoBeam2D(self, id, node_ids, E, A, ixx, G, A_s)

        fea.add_element(self, element)


class FrameElement(FiniteElement):
    """asldkjaslkd
    """

    def __init__(self, analysis, id, node_ids, E, A):
        """
        """

        super().__init__(analysis, id, node_ids)
        self.EA = E * A

    def plot_element(self, ax):
        """
        """

        coords = self.get_coords()

        ax.plot(coords[:, 0], coords[:, 1], 'k.-',
                linewidth=2, markersize=8)


class EulerBernoulliBeam2D(FrameElement):
    """
    """

    # TODO: properly implement EB with shape functions etc.

    def __init__(self, analysis, id, node_ids, E, A, ixx):
        """
        """

        super().__init__(analysis, id, node_ids, E, A)
        self.EI = E * ixx

    def get_stiff_matrix(self):
        """
        """

        coords = self.get_coords()  # coordinates of nodes

        if self.analysis.analysis_type == 'linear':
            # compute geometric parameters
            dx = coords[1]-coords[0]
            l0 = np.linalg.norm(dx)
            phi = np.arctan2(dx[1], dx[0])

            # use analytical integration
            c = np.cos(phi)
            s = np.sin(phi)

            k_el_bar = self.EA / l0 * np.array([[1, 0, 0, -1, 0, 0],
                                                [0, 0, 0, 0, 0, 0],
                                                [0, 0, 0, 0, 0, 0],
                                                [-1, 0, 0, 1, 0, 0],
                                                [0, 0, 0, 0, 0, 0],
                                                [0, 0, 0, 0, 0, 0]])

            k_el_beam = self.EI / (l0 * l0 * l0) * (
                np.array([[0, 0, 0, 0, 0, 0],
                          [0, 12, 6 * l0, 0, -12, 6 * l0],
                          [0, 6 * l0, 4 * l0 * l0, 0, -6 * l0, 2 * l0 * l0],
                          [0, 0, 0, 0, 0, 0],
                          [0, -12, -6 * l0, 0, 12, -6 * l0],
                          [0, 6 * l0, 2 * l0 * l0, 0, -6 * l0, 4 * l0 * l0]]))

            k_el = k_el_bar + k_el_beam

            T = np.array([[c, s, 0, 0, 0, 0],
                          [-s, c, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0],
                          [0, 0, 0, c, s, 0],
                          [0, 0, 0, -s, c, 0],
                          [0, 0, 0, 0, 0, 1]])

            return np.matmul(np.matmul(np.transpose(T), k_el), T)


class TimoshenkoBeam2D(FrameElement):
    """
    """

    def __init__(self, analysis, id, node_ids, E, A, ixx, G, A_s):
        """
        """

        super().__init__(analysis, id, node_ids, E, A)
        self.EI = E * ixx
        self.GA_s = G * A_s

    def get_stiff_matrix(self):
        """
        """

        coords = self.get_coords()  # coordinates of nodes

        if self.analysis.analysis_type == 'linear':
            # compute geometric parameters
            dx = coords[1]-coords[0]
            l0 = np.linalg.norm(dx)
            phi = np.arctan2(dx[1], dx[0])

            # use one point integration
            N = np.array([0.5, 0.5])

            c = np.cos(phi) / l0
            s = np.sin(phi) / l0
            t = 1 / l0

            bmat = np.array([[-c, -s, 0,  c,  s, 0],
                             [s, -c, -N[0], -s,  c, -N[1]],
                             [0,  0, -t,  0,  0, t]])

            dmat = np.diag([self.EA, self.GA_s, self.EI])

            return l0 * np.matmul(np.matmul(np.transpose(bmat), dmat), bmat)
