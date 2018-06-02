import numpy as np
from matplotlib.patches import Polygon
from fea.fea import fea, FiniteElement
from post.results import FrameForceVector
from fea.exceptions import FEAInputError


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
                    G=0, A_s=0, rho=1):
        """
        EB2: 2 Noded Euler Bernoulli Beam
        TB2: 2 Noded Timoshenko Beam
        """

        # TODO: check value types e.g. id and node_ids are ints

        if el_type == 'EB2':
            element = EulerBernoulliFrame2D(self, id, node_ids, E, A, ixx, rho)
        elif el_type == 'TB2':
            element = TimoshenkoFrame2D(self, id, node_ids, E, A, ixx, G, A_s,
                                        rho)

        fea.add_element(self, element)


class FrameElement(FiniteElement):
    """asldkjaslkd
    """

    def __init__(self, analysis, id, node_ids, E, A, rho):
        """
        """

        super().__init__(analysis, id, node_ids)
        self.EA = E * A
        self.rhoA = rho * A

    def get_geometric_properties(self):
        """
        """

        node_coords = self.get_node_coords()
        dx = node_coords[1] - node_coords[0]
        l0 = np.linalg.norm(dx)
        phi = np.arctan2(dx[1], dx[0])
        c = np.cos(phi)
        s = np.sin(phi)

        return (node_coords, dx, l0, phi, c, s)

    def get_fint(self, case_id):
        """
        """

        return self.f_int.get_result(case_id)

    def set_fint(self, case_id, f_int):
        """
        """

        (_, _, _, _, c, s) = self.get_geometric_properties()

        f = [f_int[0] * c + f_int[1] * s, -f_int[0] * s + f_int[1] * c,
             f_int[2], f_int[3] * c + f_int[4] * s,
             -f_int[3] * s + f_int[4] * c, f_int[5]]

        self.f_int.set_result(FrameForceVector(case_id, f))

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
        # TODO: IMPLEMENT DEFAULT LINEAR BEHAVIOUR

    def plot_axial_force(self, ax, case_id, scalef):
        """alskdjaklsd

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotNLine.m.
        """

        # get geometric properties
        (node_coords, dx, l0, _, _, _) = self.get_geometric_properties()

        # get internal force vector
        f_int = self.get_fint(case_id)
        n1 = -f_int.N1  # axial force at node 1 (tension positive)
        n2 = f_int.N2  # axial force at node 2 (tension positive)

        # location of node 1 and node 2
        p1 = node_coords[0, :]
        p2 = node_coords[1, :]

        # location of the axial force diagram end points
        v = np.matmul(np.array([[0, -1], [1, 0]]), dx) / l0  # direction vector
        p3 = p2 + v * scalef * n2
        p4 = p1 + v * scalef * n1

        # plot axial force line and patch
        ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0.7, 0, 0))
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0.7, 0, 0))
        ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0.7, 0, 0))
        ax.add_patch(Polygon(np.array([p1, p2, p3, p4]),
                             facecolor=(1, 0, 0), linestyle='None', alpha=0.3))

        # plot text value of axial force
        mid1 = (p1 + p4) / 2
        mid2 = (p2 + p3) / 2
        ax.text(mid1[0], mid1[1], "{:5.3g}".format(n1), size=8,
                verticalalignment='bottom')
        ax.text(mid2[0], mid2[1], "{:5.3g}".format(n2), size=8,
                verticalalignment='bottom')

    def plot_shear_force(self, ax, case_id, scalef):
        """alskdjaklsd

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotVLine.m.
        """

        # get geometric properties
        (node_coords, dx, l0, _, _, _) = self.get_geometric_properties()

        # get internal force vector
        f_int = self.get_fint(case_id)
        v1 = f_int.V1  # shear force at node 1 (cw positive)
        v2 = -f_int.V2  # shear force at node 2 (cw positive)

        # location of node 1 and node 2
        p1 = node_coords[0, :]
        p2 = node_coords[1, :]

        # location of the shear force diagram end points
        v = np.matmul(np.array([[0, -1], [1, 0]]), dx) / l0  # direction vector
        p3 = p2 + v * scalef * v2
        p4 = p1 + v * scalef * v1

        # plot shear force line and patch
        ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0, 0.3, 0))
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0, 0.3, 0))
        ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0, 0.3, 0))
        ax.add_patch(Polygon(np.array([p1, p2, p3, p4]),
                             facecolor=(0, 0.5, 0), linestyle='None',
                             alpha=0.3))

        # plot text value of shear force
        mid1 = (p1 + p4) / 2
        mid2 = (p2 + p3) / 2
        ax.text(mid1[0], mid1[1], "{:5.3g}".format(v1), size=8,
                verticalalignment='bottom')
        ax.text(mid2[0], mid2[1], "{:5.3g}".format(v2), size=8,
                verticalalignment='bottom')

    def plot_bending_moment(self, ax, case_id, scalef):
        """alskdjaklsd

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotMLine.m.
        """

        # get geometric properties
        (node_coords, dx, l0, _, _, _) = self.get_geometric_properties()

        # get internal force vector
        f_int = self.get_fint(case_id)
        m1 = f_int.M1  # bending moment at node 1 (tension positive)
        m2 = -f_int.M2  # bending moment at node 2 (tension positive)

        # location of node 1 and node 2
        p1 = node_coords[0, :]
        p2 = node_coords[1, :]

        # location of the bending moment diagram end points
        v = np.matmul(np.array([[0, -1], [1, 0]]), dx) / l0  # direction vector
        p3 = p2 + v * scalef * m2
        p4 = p1 + v * scalef * m1

        # plot bending moment line and patch
        ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0, 0, 0.7))
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0, 0, 0.7))
        ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0, 0, 0.7))
        ax.add_patch(Polygon(np.array([p1, p2, p3, p4]),
                             facecolor=(0.2, 0.4, 0.8), linestyle='None',
                             alpha=0.3))

        # plot text value of bending moment
        mid1 = (p1 + p4) / 2
        mid2 = (p2 + p3) / 2
        ax.text(mid1[0], mid1[1], "{:5.3g}".format(m1), size=8,
                verticalalignment='bottom')
        ax.text(mid2[0], mid2[1], "{:5.3g}".format(m2), size=8,
                verticalalignment='bottom')


class EulerBernoulliFrame2D(FrameElement):
    """
    """

    # TODO: properly implement EB with shape functions etc.

    def __init__(self, analysis, id, node_ids, E, A, ixx, rho):
        """
        """

        super().__init__(analysis, id, node_ids, E, A, rho)
        self.EI = E * ixx

    def get_stiff_matrix(self):
        """
        """

        if not self.analysis.non_linear:
            # compute geometric parameters
            (_, _, l0, _, c, s) = self.get_geometric_properties()

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

    def get_geometric_stiff_matrix(self, case_id):
        """
        """

        # compute geometric parameters
        (_, _, l0, _, c, s) = self.get_geometric_properties()

        # get axial force
        try:
            f_int = self.get_fint(case_id)
        except FEAInputError as error:
            print(error)

        # get axial force in element (take average of nodal values)
        N = np.mean([-f_int.N1, f_int.N2])

        # form geometric stiffness matrix
        k_el_g = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 1.2, l0/10, 0, -1.2, l0/10],
                           [0, l0/10, 2*l0*l0/15, 0, -l0/10, -l0*l0/30],
                           [0, 0, 0, 0, 0, 0],
                           [0, -1.2, -l0/10, 0, 1.2, -l0/10],
                           [0, l0/10, -l0*l0/30, 0, -l0/10, 2*l0*l0/15]])
        k_el_g *= N / l0

        # construct rotation matrix
        T = np.array([[c, s, 0, 0, 0, 0],
                      [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, -s, c, 0],
                      [0, 0, 0, 0, 0, 1]])

        return np.matmul(np.matmul(np.transpose(T), k_el_g), T)

    def get_mass_matrix(self):
        """
        """

        # compute geometric parameters
        (_, _, l0, _, c, s) = self.get_geometric_properties()

        # compute element mass matrix
        m_el = np.array([[140, 0, 0, 70, 0, 0],
                         [0, 156, 22*l0, 0, 54, -13*l0],
                         [0, 22*l0, 4*l0*l0, 0, 13*l0, -3*l0*l0],
                         [70, 0, 0, 140, 0, 0],
                         [0, 54, 13*l0, 0, 156, -22*l0],
                         [0, -13*l0, -3*l0*l0, 0, -22*l0, 4*l0*l0]])
        m_el *= self.rhoA * l0 / 420

        # construct rotation matrix
        T = np.array([[c, s, 0, 0, 0, 0],
                      [-s, c, 0, 0, 0, 0],
                      [0, 0, 1, 0, 0, 0],
                      [0, 0, 0, c, s, 0],
                      [0, 0, 0, -s, c, 0],
                      [0, 0, 0, 0, 0, 1]])

        return np.matmul(np.matmul(np.transpose(T), m_el), T)

    def plot_deformed_element(self, ax, case_id, n, def_scale, u_el=None):
        """
        """

        # nodal displacements in xy
        if u_el is None:
            u_el = self.get_nodal_displacements(case_id)

        # compute frame geometric parameters
        (node_coords, _, l0, _, c, s) = self.get_geometric_properties()

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
            u_x[i] = u * c - v * s
            u_y[i] = u * s + v * c

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

    def __init__(self, analysis, id, node_ids, E, A, ixx, G, A_s, rho):
        """
        """

        super().__init__(analysis, id, node_ids, E, A, rho)
        self.EI = E * ixx
        self.GA_s = G * A_s
