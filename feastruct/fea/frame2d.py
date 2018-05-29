import numpy as np
from fea.fea import fea, FiniteElement


class Frame2D(fea):
    """
    """

    def __init__(self, nodes=None, elements=None, freedom_cases=None,
                 load_cases=None, analysis_cases=None, non_linear=False):
        """
        """

        # initialise parent fea class
        super().__init__(nodes, elements, freedom_cases, load_cases,
                         analysis_cases, non_linear, dofs=3)

    def add_element(self, id, node_ids, el_type='EB2', E=1, A=1, ixx=1,
                    G=0, A_s=0):
        """
        EB2: 2 Noded Euler Bernoulli Beam
        TB2: 2 Noded Timoshenko Beam
        """

        # TODO: check value types e.g. id and node_ids are ints

        if el_type == 'EB2':
            element = EulerBernoulliFrame2D(self, id, node_ids, E, A, ixx)
        elif el_type == 'TB2':
            element = TimoshenkoFrame2D(self, id, node_ids, E, A, ixx, G, A_s)

        fea.add_element(self, element)


class FrameElement(FiniteElement):
    """asldkjaslkd
    """

    def __init__(self, analysis, id, node_ids, E, A):
        """
        """

        super().__init__(analysis, id, node_ids)
        self.EA = E * A

    def plot_element(self, ax, linestyle='-', linewidth=2, marker='.'):
        """
        """

        coords = self.get_node_coords()

        ax.plot(coords[:, 0], coords[:, 1], 'k.', linestyle=linestyle,
                linewidth=linewidth, marker=marker, markersize=8)

    def plot_deformed_element(self, ax):
        """
        """
        pass
        # IMPLEMENT DEFAULT LINEAR BEHAVIOUR


class EulerBernoulliFrame2D(FrameElement):
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

        node_coords = self.get_node_coords()  # coordinates of nodes

        if not self.analysis.non_linear:
            # compute geometric parameters
            dx = node_coords[1]-node_coords[0]
            l0 = np.linalg.norm(dx)
            phi = np.arctan2(dx[1], dx[0])
            c = np.cos(phi)
            s = np.sin(phi)

            # use analytical integration result:
            # compute bar stiffness
            k_el_bar = self.EA / l0 * np.array([[1, 0, 0, -1, 0, 0],
                                                [0, 0, 0, 0, 0, 0],
                                                [0, 0, 0, 0, 0, 0],
                                                [-1, 0, 0, 1, 0, 0],
                                                [0, 0, 0, 0, 0, 0],
                                                [0, 0, 0, 0, 0, 0]])

            # compute beam stiffness
            k_el_beam = self.EI / (l0 * l0 * l0) * (
                np.array([[0, 0, 0, 0, 0, 0],
                          [0, 12, 6 * l0, 0, -12, 6 * l0],
                          [0, 6 * l0, 4 * l0 * l0, 0, -6 * l0, 2 * l0 * l0],
                          [0, 0, 0, 0, 0, 0],
                          [0, -12, -6 * l0, 0, 12, -6 * l0],
                          [0, 6 * l0, 2 * l0 * l0, 0, -6 * l0, 4 * l0 * l0]]))

            k_el = k_el_bar + k_el_beam

            # construct rotation matrix
            T = np.array([[c, s, 0, 0, 0, 0],
                          [-s, c, 0, 0, 0, 0],
                          [0, 0, 1, 0, 0, 0],
                          [0, 0, 0, c, s, 0],
                          [0, 0, 0, -s, c, 0],
                          [0, 0, 0, 0, 0, 1]])

            return np.matmul(np.matmul(np.transpose(T), k_el), T)

    def plot_deformed_element(self, ax, case_id, n, def_scale):
        """
        """

        node_coords = self.get_node_coords()  # coordinates of nodes
        # nodal displacements in xy
        u_el = self.get_nodal_displacements(case_id)

        # compute frame geometric parameters
        dx = node_coords[1]-node_coords[0]
        l0 = np.linalg.norm(dx)
        phi = np.arctan2(dx[1], dx[0])
        c = np.cos(phi)
        s = np.sin(phi)

        # rotate nodal displacements to local axis
        T = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
        u_el[0, :] = np.matmul(T, u_el[0, :])
        u_el[1, :] = np.matmul(T, u_el[1, :])

        # allocate displacement vectors
        u_x = np.zeros(n)
        u_y = np.zeros(n)
        x = np.zeros(n)
        y = np.zeros(n)

        # original location of frame station points
        x0 = np.linspace(node_coords[0, 0], node_coords[1, 0], n)
        y0 = np.linspace(node_coords[0, 1], node_coords[1, 1], n)

        # loop through stations on frame
        for (i, xi) in enumerate(np.linspace(-1, 1, n)):
            # element shape functions
            N_u = np.array([0.5 - xi / 2, 0.5 + xi / 2])
            N_v = np.array([0.25 * (1 - xi) * (1 - xi) * (2 + xi),
                            0.125 * l0 * (1 - xi) * (1 - xi) * (1 + xi),
                            0.25 * (1 + xi) * (1 + xi) * (2 - xi),
                            -0.125 * l0 * (1 + xi) * (1 + xi) * (1 - xi)])

            # compute local displacements
            u = np.dot(N_u, np.array([u_el[0, 0], u_el[1, 0]]))
            v = np.dot(N_v, np.array([u_el[0, 1], u_el[0, 2],
                                      u_el[1, 1], u_el[1, 2]]))

            # scale displacements by deformation scale
            u *= def_scale
            v *= def_scale

            # compute cartesian displacements at point i
            u_x[i] = u * np.cos(phi) - v * np.sin(phi)
            u_y[i] = u * np.sin(phi) + v * np.cos(phi)

            # compute location of point i
            x[i] = u_x[i] + x0[i]
            y[i] = u_y[i] + y0[i]

        # plot frame elements
        for i in range(n - 1):
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], 'k-', linewidth=2)

        # plot end markers
        ax.plot(x[0], y[0], 'k.', markersize=8)
        ax.plot(x[-1], y[-1], 'k.', markersize=8)


class TimoshenkoFrame2D(FrameElement):
    """ TODO: implement with proper shape functions
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

        coords = self.get_node_coords()  # coordinates of nodes

        if not self.analysis.non_linear:
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
