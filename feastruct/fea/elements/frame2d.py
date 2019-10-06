import numpy as np
from scipy import integrate
from matplotlib.patches import Polygon
from feastruct.fea.elements.frame import FrameElement
from feastruct.fea.bcs import ElementLoad


class FrameElement2D(FrameElement):
    """Class for a 2D frame element.

    Provides a number of methods that can be used for a 2D frame element.

    :cvar nodes: List of node objects defining the element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
    :cvar efs: Element freedom signature
    :vartype efs: list[bool]
    :cvar f_int: List of internal force vector results stored for each analysis case
    :vartype f_int: list[:class:`~feastruct.fea.fea.ForceVector`]
    :cvar section: Section object for the element
    :vartype section: :class:`~feastruct.pre.section.Section`
    """

    def __init(self, nodes, material, efs, section):
        """Inits the FrameElement2D class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :cvar efs: Element freedom signature
        :vartype efs: list[bool]
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # initialise parent FrameElement class
        super().__init__(nodes=nodes, material=material, efs=efs, section=section)

    def plot_element(self, ax, linestyle='-', linewidth=2, marker='.'):
        """Plots the undeformed frame element on the axis defined by ax.

        :param ax: Axis object on which to draw the element
        :type ax: :class:`matplotlib.axes.Axes`
        :param string linestyle: Element linestyle
        :param int linewidth: Element linewidth
        :param string marker: Node marker type
        """

        coords = self.get_node_coords()

        ax.plot(coords[:, 0], coords[:, 1], 'k.', linestyle=linestyle,
                linewidth=linewidth, marker=marker, markersize=8)

    def plot_deformed_element(self, ax, analysis_case, n, def_scale, u_el=None):
        """Plots a 2D frame element in its deformed configuration for the displacement vector
        defined by the analysis_case. The deformation is based on the element shape functions. If a
        displacement vector, *u_el*, is supplied, uses this to plot the deformed element.

        :param ax: Axis object on which to draw the element
        :type ax: :class:`matplotlib.axes.Axes`
        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int n: Number of linear subdivisions used to plot the element
        :param float def_scale: Deformation scale
        :param u_el: Element displacement vector of size *(n_node x n_dof)*
        :type u_el: :class:`numpy.ndarray`
        """

        # if no displacement vector is supplied, get the analysis results
        if u_el is None:
            # get local axis displacements
            disps = self.get_displacements(n=n, analysis_case=analysis_case)
            xis = disps[:, 0]
            us = disps[:, 1]
            vs = disps[:, 2]

        # if a displacement vector is supplied
        else:
            # set stations
            xis = np.linspace(-1, 1, n)

            # rotate nodal displacements to local axis
            T = self.get_transformation_matrix()
            u_el[0, :] = np.matmul(T, u_el[0, :])
            u_el[1, :] = np.matmul(T, u_el[1, :])

        # redefine number of stations
        n = len(xis)

        # compute frame geometric parameters
        (node_coords, _, _, c) = self.get_geometric_properties()

        # allocate displacement vectors
        u_x = np.zeros(n)
        u_y = np.zeros(n)
        x = np.zeros(n)
        y = np.zeros(n)

        # original location of frame station points
        x0 = np.linspace(node_coords[0, 0], node_coords[1, 0], n)
        y0 = np.linspace(node_coords[0, 1], node_coords[1, 1], n)

        # loop through stations on frame
        for (i, xi) in enumerate(xis):
            if u_el is None:
                u = us[i]
                v = vs[i]
            else:
                (u, v, _) = self.calculate_local_displacement(xi, u_el)

            # scale displacements by deformation scale
            u *= def_scale
            v *= def_scale

            # compute cartesian displacements at point i
            u_x[i] = u * c[0] - v * c[1]
            u_y[i] = u * c[1] + v * c[0]

            # compute location of point i
            x[i] = u_x[i] + x0[i]
            y[i] = u_y[i] + y0[i]

        # plot frame elements
        for i in range(n - 1):
            ax.plot([x[i], x[i+1]], [y[i], y[i+1]], 'k-', linewidth=2)

        # plot end markers
        ax.plot(x[0], y[0], 'k.', markersize=8)
        ax.plot(x[-1], y[-1], 'k.', markersize=8)

    def plot_axial_force(self, ax, analysis_case, scalef, n):
        """Plots the axial force diagram from a static analysis defined by case_id. N.B. this
        method is adapted from the MATLAB code by F.P. van der Meer: plotNLine.m.

        :param ax: Axis object on which to draw the element
        :type ax: :class:`matplotlib.axes.Axes`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param float scalef: Factor by which to scale the axial force diagram
        :param int n: Number of points at which to plot the axial force diagram
        """

        # get geometric properties
        (node_coords, dx, l0, _) = self.get_geometric_properties()

        # get axial force diagram
        (xis, afd) = self.get_afd(n=n, analysis_case=analysis_case)

        # get indices of min and max values of bending moment
        min_index = np.argmin(afd)
        max_index = np.argmax(afd)

        # get end node coordinates
        end1 = node_coords[0, 0:2]
        end2 = node_coords[1, 0:2]

        # plot shear force diagram
        for (i, xi) in enumerate(xis[:-1]):
            n1 = afd[i]
            n2 = afd[i+1]

            a_i1 = (xi + 1) / 2  # percentage along element for xis[i]
            a_i2 = (xis[i+1] + 1) / 2  # percentage along element for xis[i+1]

            # location of node 1 and node 2
            p1 = end1 + a_i1 * (end2 - end1)
            p2 = end1 + a_i2 * (end2 - end1)

            # location of the axial force diagram end points
            v = np.matmul(np.array([[0, -1], [1, 0]]), dx[0:2]) / l0  # direction vector
            p3 = p2 + v * scalef * n2
            p4 = p1 + v * scalef * n1

            # plot shear force line and patch
            ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0.7, 0, 0))
            ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0.7, 0, 0))
            ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0.7, 0, 0))
            ax.add_patch(Polygon(
                np.array([p1, p2, p3, p4]), facecolor=(1, 0, 0), linestyle='None', alpha=0.3
            ))

            # plot end text values of bending moment
            if i == 0:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(n1), size=8, verticalalignment='bottom')
            elif i == len(xis) - 2:
                mid = (p2 + p3) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(n2), size=8, verticalalignment='bottom')

            # plot text value of min bending moment
            if i == min_index:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(n1), size=8, verticalalignment='bottom')

            # plot text value of max bending moment
            if i == max_index:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(n1), size=8, verticalalignment='bottom')

    def plot_shear_force(self, ax, analysis_case, scalef, n):
        """Plots the axial force diagram from a static analysis defined by case_id. N.B. this
        method is adapted from the MATLAB code by F.P. van der Meer: plotVLine.m.

        :param ax: Axis object on which to draw the element
        :type ax: :class:`matplotlib.axes.Axes`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param float scalef: Factor by which to scale the shear force diagram
        :param int n: Number of points at which to plot the shear force diagram
        """

        # get geometric properties
        (node_coords, dx, l0, _) = self.get_geometric_properties()

        # get shear force diagram
        (xis, sfd) = self.get_sfd(n=n, analysis_case=analysis_case)

        # get indices of min and max values of bending moment
        min_index = np.argmin(sfd)
        max_index = np.argmax(sfd)

        # get end node coordinates
        end1 = node_coords[0, 0:2]
        end2 = node_coords[1, 0:2]

        # plot shear force diagram
        for (i, xi) in enumerate(xis[:-1]):
            v1 = sfd[i]
            v2 = sfd[i+1]

            a_i1 = (xi + 1) / 2  # percentage along element for xis[i]
            a_i2 = (xis[i+1] + 1) / 2  # percentage along element for xis[i+1]

            # location of node 1 and node 2
            p1 = end1 + a_i1 * (end2 - end1)
            p2 = end1 + a_i2 * (end2 - end1)

            # location of the shear force diagram end points
            v = np.matmul(np.array([[0, -1], [1, 0]]), dx[0:2]) / l0  # direction vector
            p3 = p2 + v * scalef * v2
            p4 = p1 + v * scalef * v1

            # plot shear force line and patch
            ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0, 0.3, 0))
            ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0, 0.3, 0))
            ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0, 0.3, 0))
            ax.add_patch(Polygon(
                np.array([p1, p2, p3, p4]), facecolor=(0, 0.5, 0), linestyle='None', alpha=0.3
            ))

            # plot end text values of bending moment
            if i == 0:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(v1), size=8, verticalalignment='bottom')
            elif i == len(xis) - 2:
                mid = (p2 + p3) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(v2), size=8, verticalalignment='bottom')

            # plot text value of min bending moment
            if i == min_index:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(v1), size=8, verticalalignment='bottom')

            # plot text value of max bending moment
            if i == max_index:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(v1), size=8, verticalalignment='bottom')

    def plot_bending_moment(self, ax, analysis_case, scalef, n):
        """Plots the axial force diagram from a static analysis defined by case_id. N.B. this
        method is adapted from the MATLAB code by F.P. van der Meer: plotMLine.m.

        :param ax: Axis object on which to draw the element
        :type ax: :class:`matplotlib.axes.Axes`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param float scalef: Factor by which to scale the bending moment diagram
        :param int n: Number of points at which to plot the bending moment diagram
        """

        # get geometric properties
        (node_coords, dx, l0, _) = self.get_geometric_properties()

        # get bending moment diagram
        (xis, bmd) = self.get_bmd(n=n, analysis_case=analysis_case)

        # get indices of min and max values of bending moment
        min_index = np.argmin(bmd)
        max_index = np.argmax(bmd)

        # get end node coordinates
        end1 = node_coords[0, 0:2]
        end2 = node_coords[1, 0:2]

        # plot bending moment diagram
        for (i, xi) in enumerate(xis[:-1]):
            m1 = bmd[i]
            m2 = bmd[i+1]

            a_i1 = (xi + 1) / 2  # percentage along element for xis[i]
            a_i2 = (xis[i+1] + 1) / 2  # percentage along element for xis[i+1]

            # location of node 1 and node 2
            p1 = end1 + a_i1 * (end2 - end1)
            p2 = end1 + a_i2 * (end2 - end1)

            # location of the bending moment diagram end points
            v = np.matmul(np.array([[0, -1], [1, 0]]), dx[0:2]) / l0  # direction vector
            p3 = p2 + v * scalef * m2
            p4 = p1 + v * scalef * m1

            # plot bending moment line and patch
            ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0, 0, 0.7))
            ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0, 0, 0.7))
            ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0, 0, 0.7))
            ax.add_patch(Polygon(
                np.array([p1, p2, p3, p4]), facecolor=(0.2, 0.4, 0.8), linestyle='None', alpha=0.3
            ))

            # plot end text values of bending moment
            if i == 0:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(m1), size=8, verticalalignment='bottom')
            elif i == len(xis) - 2:
                mid = (p2 + p3) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(m2), size=8, verticalalignment='bottom')

            # plot text value of min bending moment
            if i == min_index:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(m1), size=8, verticalalignment='bottom')

            # plot text value of max bending moment
            if i == max_index:
                mid = (p1 + p4) / 2
                ax.text(mid[0], mid[1], "{:.3e}".format(m1), size=8, verticalalignment='bottom')


class Bar2D_2N(FrameElement2D):
    """Two noded, two dimensional bar element that can resist an axial force only. The element is
    defined by its two end nodes and uses two linear shape functions to obtain analytical results.

    :cvar nodes: List of node objects defining the element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
    :cvar efs: Element freedom signature
    :vartype efs: list[bool]
    :cvar f_int: List of internal force vector results stored for each analysis case
    :vartype f_int: list[:class:`~feastruct.fea.fea.ForceVector`]
    :cvar section: Section object for the element
    :vartype section: :class:`~feastruct.pre.section.Section`
    """

    def __init__(self, nodes, material, section):
        """Inits the Bar2D_2N class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # set the element freedom signature
        efs = [True, True, False, False, False, False]

        # initialise parent FrameElement2D class
        super().__init__(nodes=nodes, material=material, efs=efs, section=section)

    def get_shape_function(self, xi):
        """Returns the value of the shape functions *N1* and *N2* at *xi*.

        :param float xi: Position along the element

        :returns: Value of the shape functions *(N1, N2)* at *xi*
        :rtype: :class:`numpy.ndarray`
        """

        return np.array([0.5 - xi / 2, 0.5 + xi / 2])

    def get_stiffness_matrix(self):
        """Gets the stiffness matrix for a two noded, 2D bar element. The stiffness matrix has been
        analytically integrated so numerical integration is not necessary.

        :returns: 4 x 4 element stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        E = self.material.elastic_modulus
        A = self.section.area
        cx = c[0]
        cy = c[1]

        # construct rotation matrix
        T = np.array([
            [cx, cy, 0, 0],
            [0, 0, cx, cy]
        ])

        # compute bar stiffness matrix
        k = E * A / l0 * np.array([
            [1, -1],
            [-1, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), k), T)

    def get_geometric_stiff_matrix(self, analysis_case):
        """Gets the geometric stiffness matrix for a two noded, 2D bar element. The stiffness
        matrix has been analytically integrated so numerical integration is not necessary. The
        geometric stiffness matrix requires an axial force so the analysis_case from a static
        analysis must be provided.

        :param analysis_case: Analysis case from which to extract the axial force
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: 4 x 4 element geometric stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        cx = c[0]
        cy = c[1]

        # get axial force
        f_int = self.get_fint(analysis_case)

        # get axial force in element (take average of nodal values)
        N = np.mean([-f_int[0], f_int[1]])

        # construct rotation matrix
        T = np.array([
            [0, 0, cx, cy],
            [cx, cy, 0, 0]
        ])

        # compute bar geometric stiffness matrix
        k_g = N / l0 * np.array([
            [1, -1],
            [-1, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), k_g), T)

    def get_mass_matrix(self):
        """Gets the mass matrix for a for a two noded, 2D bar element. The mass matrix has been
        analytically integrated so numerical integration is not necessary.

        :returns: 4 x 4 element mass matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        rho = self.material.rho
        A = self.section.area
        cx = c[0]
        cy = c[1]

        T = np.array([
            [cx, cy, 0, 0],
            [0, 0, cx, cy]
        ])

        # compute element mass matrix
        m_el = rho * A * l0 / 6 * np.array([
            [2, 1],
            [1, 2]
        ])

        return np.matmul(np.matmul(np.transpose(T), m_el), T)

    def get_internal_actions(self, analysis_case):
        """Returns the internal actions for a two noded 2D bar element.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: An array containing the internal actions for the element
            *(N1, N2)*
        :rtype: :class:`numpy.ndarray`
        """

        (_, _, _, c) = self.get_geometric_properties()
        f_int = self.get_fint(analysis_case=analysis_case)

        cx = c[0]
        cy = c[1]

        f = np.array([
            f_int[0] * cx + f_int[1] * cy,
            f_int[2] * cx + f_int[3] * cy
        ])

        return f

    def get_displacements(self, n, analysis_case):
        """Returns a list of the local displacements, *(u, v, w, ru, rv, rw)*, along the element
        for the analysis case and a minimum of *n* subdivisions. A list of the stations, *xi*, is
        also included.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int n: Minimum number of sampling points

        :returns: 2D numpy array containing stations and local displacements. Each station contains
            an array of the following format: *[xi, u, v, w, ru, rv, rw]*
        :rtype: :class:`numpy.ndarray`
        """

        # get a list of the stations
        stations = self.get_sampling_points(n=n, analysis_case=analysis_case)

        # allocate results vector
        results = np.zeros((len(stations), 7))

        # loop through all the stations
        for (i, xi) in enumerate(stations):
            # get the nodal displacements
            u_el = self.get_nodal_displacements(analysis_case)

            # rotate nodal displacements to local axis
            T = self.get_transformation_matrix()
            u_el_local = np.transpose(np.matmul(T, np.transpose(u_el)))

            # element shape function at station location
            N = self.get_shape_function(xi)

            # compute local displacements
            u = np.dot(N, np.array([u_el_local[0, 0], u_el_local[1, 0]]))
            v = np.dot(N, np.array([u_el_local[0, 1], u_el_local[1, 1]]))

            # save results
            results[i, 0] = xi  # station location
            results[i, 1] = u  # u-translation
            results[i, 2] = v  # v-translation
            results[i, 3:] = None  # other dofs are not assigned

        return results

    def get_transformation_matrix(self):
        """Returns the transformation matrix for a Bar2D_2N element.

        :returns: Element transformation matrix
        :rtype: :class:`numpy.ndarray`
        """

        (_, _, _, c) = self.get_geometric_properties()

        return np.array([
            [c[0], c[1]],
            [-c[1], c[0]]
        ])

    def get_afd(self, n, analysis_case):
        """Returns the axial force diagram within the element for *n* stations for an
        analysis_case.

        :param int n: Number of stations to sample the axial force diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and axial force diagram, *afd* - *(xis, afd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # get internal forces
        f = self.get_internal_actions(analysis_case=analysis_case)
        N1 = -f[0]
        N2 = f[1]

        # allocate the axial force diagram
        afd = np.zeros(n)

        # generate list of stations
        stations = self.get_sampling_points(n=n, analysis_case=analysis_case)

        # loop over each station
        for (i, xi) in enumerate(stations):
            # get shape functions at xi
            N = self.get_shape_function(xi)

            # compute local displacements
            afd[i] = np.dot(N, np.array([N1, N2]))

        return (stations, afd)

    def get_sfd(self, n, analysis_case):
        """Returns the shear force diagram within the element for *n* stations for an
        analysis_case.

        :param int n: Number of stations to sample the shear force diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and shear force diagram, *sfd* - *(xis, sfd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # no shear force in this element
        return (np.linspace(-1, 1, n), np.zeros(n))

    def get_bmd(self, n, analysis_case):
        """Returns the bending moment diagram within the element for *n* stations for an
        analysis_case.

        :param int n: Number of stations to sample the bending moment diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and bending moment diagram, *bmd* - *(xis, bmd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # no bending moment in this element
        return (np.linspace(-1, 1, n), np.zeros(n))

    def calculate_local_displacement(self, xi, u_el):
        """Calculates the local displacement of the element at position *xi* given the displacement
        vector *u_el*.

        :param float xi: Position along the element
        :param u_el: Element displacement vector
        :type u_el: :class:`numpy.ndarray`

        :returns: Local displacement of the element *(u, v, w)*
        :rtype: tuple(float, float, float)
        """

        # element shape function
        N = self.get_shape_function(xi)

        # compute local displacements
        u = np.dot(N, np.array([u_el[0, 0], u_el[1, 0]]))
        v = np.dot(N, np.array([u_el[0, 1], u_el[1, 1]]))

        return (u, v, None)


class EulerBernoulli2D_2N(FrameElement2D):
    """Two noded, two dimensional frame element based on the Euler-Bernoulli beam formulation for
    relatively thin beams. The element is defined by its two end nodes and uses four cubic
    polynomial shape functions to obtain analytical results.

    :cvar nodes: List of node objects defining the element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
    :cvar efs: Element freedom signature
    :vartype efs: list[bool]
    :cvar f_int: List of internal force vector results stored for each analysis case
    :vartype f_int: list[:class:`~feastruct.fea.fea.ForceVector`]
    :cvar section: Section object for the element
    :vartype section: :class:`~feastruct.pre.section.Section`
    """

    def __init__(self, nodes, material, section):
        """Inits the EulerBernoulli2D_2N class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # set the element freedom signature
        efs = [True, True, False, False, False, True]

        # initialise parent FrameElement2D class
        super().__init__(nodes=nodes, material=material, efs=efs, section=section)

    def get_shape_function(self, xi):
        """Returns the value of the shape functions *Nu1*, *Nu2*, *Nv1* and *Nv2* at *xi*.

        :param float xi: Position along the element

        :returns: Value of the shape functions *((Nu1, Nu2), (Nv1, Nv2))* at *xi*
        :rtype: :class:`numpy.ndarray`
        """

        # compute frame geometric parameters
        (_, _, l0, _) = self.get_geometric_properties()

        # element shape functions
        N_u = np.array([0.5 - xi / 2, 0.5 + xi / 2])
        N_v = np.array([
            0.25 * (1 - xi) * (1 - xi) * (2 + xi),
            0.125 * l0 * (1 - xi) * (1 - xi) * (1 + xi),
            0.25 * (1 + xi) * (1 + xi) * (2 - xi),
            -0.125 * l0 * (1 + xi) * (1 + xi) * (1 - xi)
        ])

        return (N_u, N_v)

    def get_stiffness_matrix(self):
        """Gets the stiffness matrix for a two noded 2D Euler-Bernoulli frame element. The
        stiffness matrix has been analytically integrated so numerical integration is not
        necessary.

        :returns: 6 x 6 element stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        E = self.material.elastic_modulus
        A = self.section.area
        ixx = self.section.ixx
        cx = c[0]
        cy = c[1]

        # compute bar stiffness matrix
        k_el_bar = E * A / l0 * np.array([
            [1, 0, 0, -1, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [-1, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0]
        ])

        # compute beam stiffness matrix
        k_el_beam = E * ixx / (l0 * l0 * l0) * np.array([
            [0, 0, 0, 0, 0, 0],
            [0, 12, 6*l0, 0, -12, 6*l0],
            [0, 6*l0, 4*l0*l0, 0, -6*l0, 2*l0*l0],
            [0, 0, 0, 0, 0, 0],
            [0, -12, -6*l0, 0, 12, -6*l0],
            [0, 6*l0, 2*l0*l0, 0, -6*l0, 4*l0*l0]
        ])

        k_el = k_el_bar + k_el_beam

        # construct rotation matrix
        T = np.array([
            [cx, cy, 0, 0, 0, 0],
            [-cy, cx, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, cx, cy, 0],
            [0, 0, 0, -cy, cx, 0],
            [0, 0, 0, 0, 0, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), k_el), T)

    def get_geometric_stiff_matrix(self, analysis_case):
        """Gets the geometric stiffness matrix for a two noded 2D Euler-Bernoulli frame element.
        The stiffness matrix has been analytically integrated so numerical integration is not
        necessary. The geometric stiffness matrix requires an axial force so the analysis_case from
        a static analysis must be provided.

        :param analysis_case: Analysis case from which to extract the axial force
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: 6 x 6 element geometric stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        cx = c[0]
        cy = c[1]

        # get axial force
        f_int = self.get_fint(analysis_case)

        # get axial force in element (take average of nodal values)
        N = np.mean([-f_int[0], f_int[3]])

        # form geometric stiffness matrix
        k_el_g = np.array([
            [0, 0, 0, 0, 0, 0],
            [0, 1.2, 0.1*l0, 0, -1.2, 0.1*l0],
            [0, 0.1*l0, 2*l0*l0/15.0, 0, -0.1*l0, -l0*l0/30.0],
            [0, 0, 0, 0, 0, 0],
            [0, -1.2, -0.1*l0, 0, 1.2, -0.1*l0],
            [0, 0.1*l0, -l0*l0/30.0, 0, -0.1*l0, 2*l0*l0/15.0]
        ])
        k_el_g *= N / l0

        # construct rotation matrix
        T = np.array([
            [cx, cy, 0, 0, 0, 0],
            [-cy, cx, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, cx, cy, 0],
            [0, 0, 0, -cy, cx, 0],
            [0, 0, 0, 0, 0, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), k_el_g), T)

    def get_mass_matrix(self):
        """Gets the mass matrix for a for a two noded 2D Euler-Bernoulli frame element. The mass
        matrix has been analytically integrated so numerical integration is not necessary.

        :returns: 6 x 6 element mass matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        rho = self.material.rho
        A = self.section.area
        cx = c[0]
        cy = c[1]

        # compute element mass matrix
        m_el = np.array([
            [140, 0, 0, 70, 0, 0],
            [0, 156, 22*l0, 0, 54, -13*l0],
            [0, 22*l0, 4*l0*l0, 0, 13*l0, -3*l0*l0],
            [70, 0, 0, 140, 0, 0],
            [0, 54, 13*l0, 0, 156, -22*l0],
            [0, -13*l0, -3*l0*l0, 0, -22*l0, 4*l0*l0]
        ])
        m_el *= rho * A * l0 / 420

        # construct rotation matrix
        T = np.array([
            [cx, cy, 0, 0, 0, 0],
            [-cy, cx, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, cx, cy, 0],
            [0, 0, 0, -cy, cx, 0],
            [0, 0, 0, 0, 0, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), m_el), T)

    def get_internal_actions(self, analysis_case):
        """Returns the internal actions for a two noded 2D Euler-Bernoulli frame element.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: An array containing the internal actions for the element
            *(N1, V1, M1, N2, V2, M2)*
        :rtype: :class:`numpy.ndarray`
        """

        (_, _, _, c) = self.get_geometric_properties()
        f_int = self.get_fint(analysis_case=analysis_case)

        cx = c[0]
        cy = c[1]

        f = np.array([
            f_int[0] * cx + f_int[1] * cy,
            -f_int[0] * cy + f_int[1] * cx,
            f_int[2],
            f_int[3] * cx + f_int[4] * cy,
            -f_int[3] * cy + f_int[4] * cx,
            f_int[5]
        ])

        return f

    def get_displacements(self, n, analysis_case):
        """Returns a list of the local displacements, *(u, v, w, ru, rv, rw)*, along the element
        for the analysis case and a minimum of *n* subdivisions. A list of the stations, *xi*, is
        also included.

        An extra station is included if there is a point of zero rotation resulting in a local
        displacement maxima/minima.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int n: Minimum number of sampling points

        :returns: 2D numpy array containing stations and local displacements. Each station contains
            an array of the following format: *[xi, u, v, w, ru, rv, rw]*
        :rtype: :class:`numpy.ndarray`
        """

        # get a list of the stations
        stations = self.get_sampling_points(n=n, analysis_case=analysis_case, defl=True)

        # allocate results vector
        results = np.zeros((len(stations), 7))

        # loop through all the stations
        for (i, xi) in enumerate(stations):
            # get the nodal displacements
            u_el = self.get_nodal_displacements(analysis_case)

            # rotate nodal displacements to local axis
            T = self.get_transformation_matrix()
            u_el_local = np.transpose(np.matmul(T, np.transpose(u_el)))

            # element axial shape function at station location
            (N_u, _) = self.get_shape_function(xi)

            # compute axial displacement
            u = np.dot(N_u, np.array([u_el_local[0, 0], u_el_local[1, 0]]))

            # calculate transverse displacement
            v = self.calculate_transverse_displacement(
                xis=[xi], v0=u_el_local[0, 1], phi0=u_el_local[0, 2],
                analysis_case=analysis_case)[0]

            # calculate rotation
            rz = self.calculate_rotation(
                xis=[xi], phi0=u_el_local[0, 2], analysis_case=analysis_case)[0]

            # save results
            results[i, 0] = xi  # station location
            results[i, 1] = u  # u-translation
            results[i, 2] = v  # v-translation
            results[i, 5] = rz  # w-rotation
            results[i, 3:5] = None  # other dofs are not assigned

        return results

    def calculate_rotation(self, xis, phi0, analysis_case):
        """integrate bending moment to get rotations

        vectorised to allow for fixed quadrature
        """

        # ensure input is an array
        if type(xis) is float:
            xis = np.array([xis])

        # allocate rotations
        phis = np.zeros(len(xis))

        # get bending stiffness
        E = self.material.elastic_modulus
        ixx = self.section.ixx

        # curvature function: kappa(x) = -M(x) / EI
        def kappa(x): return -self.get_bm(x, analysis_case) / E / ixx

        for (i, xi) in enumerate(xis):
            # integrate curvature to get change in rotation from x = -1 to x = xi
            (delta_phi, err) = integrate.fixed_quad(kappa, -1, xi)

            # get member length
            (_, _, l0, _) = self.get_geometric_properties()

            # transform xi to l0
            delta_phi *= l0 / 2

            phis[i] = phi0 + delta_phi

        # rotation = initial rotation + change in rotation
        return phis

    def calculate_transverse_displacement(self, xis, v0, phi0, analysis_case):
        """integrate rotations to get transvere displacement

        vectorised to allow for fixed quadrature
        """

        # allocate displacements
        disps = np.zeros(len(xis))

        # rotation function
        def phi(x): return self.calculate_rotation(x, phi0, analysis_case)

        for (i, xi) in enumerate(xis):
            # integrate rotation to get change in transverse displacement from x = -1 to x = xi
            (delta_v, err) = integrate.fixed_quad(phi, -1, xi)

            # get member length
            (_, _, l0, _) = self.get_geometric_properties()

            # transform xi to l0
            delta_v *= l0 / 2

            # transvere disp = initial transvere disp + change in transvere disp
            disps[i] = v0 + delta_v

        return disps

    def get_transformation_matrix(self):
        """Returns the transformation matrix for an EulerBernoulli2D_2N element.

        :returns: Element transformation matrix
        :rtype: :class:`numpy.ndarray`
        """

        (_, _, _, c) = self.get_geometric_properties()

        return np.array([[c[0], c[1], 0], [-c[1], c[0], 0], [0, 0, 1]])

    def get_afd(self, n, analysis_case):
        """Returns the axial force diagram within the element for a minimum of *n* stations for an
        analysis_case.

        :param int n: Minimum number of stations to sample the axial force diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and axial force diagram, *afd* - *(xis, afd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # get internal forces
        f = self.get_internal_actions(analysis_case=analysis_case)
        N1 = -f[0]
        N2 = f[3]

        # generate list of stations
        stations = self.get_sampling_points(n=n, analysis_case=analysis_case)

        # allocate the axial force diagram
        afd = np.zeros(len(stations))

        # loop over each station
        for (i, xi) in enumerate(stations):
            # get shape functions at xi
            (N, _) = self.get_shape_function(xi)

            # compute axial force diagram
            afd[i] = np.dot(N, np.array([N1, N2]))

        return (stations, afd)

    def get_sfd(self, n, analysis_case):
        """Returns the shear force diagram within the element for a minimum of *n* stations for an
        analysis_case.

        :param int n: Minimum number of stations to sample the shear force diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and shear force diagram, *sfd* - *(xis, sfd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # get internal forces
        f = self.get_internal_actions(analysis_case=analysis_case)
        V1 = f[1]
        V2 = -f[4]

        # get sampling points
        stations = self.get_sampling_points(n=n, analysis_case=analysis_case)

        # get list of element loads
        element_loads = self.get_element_loads(analysis_case=analysis_case)

        # allocate the shear force diagram
        sfd = np.zeros(len(stations))

        # loop over each station
        for (i, xi) in enumerate(stations):
            # get shape functions at xi
            (N, _) = self.get_shape_function(xi)

            # compute shear force diagram
            sfd[i] = np.dot(N, np.array([V1, V2]))

            # add shear force due to element loads
            for element_load in element_loads:
                sfd[i] += element_load.get_internal_sfd(xi)

        return (stations, sfd)

    def get_bmd(self, n, analysis_case):
        """Returns the bending moment diagram within the element for *n* stations for an
        analysis_case. An additional station is added at all locations where the shear force is
        zero to ensure that bending moment maxima/minima are captured.

        :param int n: Number of stations to sample the bending moment diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and bending moment diagram, *bmd* - *(xis, bmd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        # get internal forces
        f = self.get_internal_actions(analysis_case=analysis_case)
        M1 = f[2]
        M2 = -f[5]

        # get sampling points
        stations = self.get_sampling_points(n=n, analysis_case=analysis_case, bm=True)

        # get list of element loads
        element_loads = self.get_element_loads(analysis_case=analysis_case)

        # allocate the bending moment diagram
        bmd = np.zeros(len(stations))

        # loop over each station
        for (i, xi) in enumerate(stations):
            # get shape functions at xi
            (N, _) = self.get_shape_function(xi)

            # compute bending moment diagram
            bmd[i] = np.dot(N, np.array([M1, M2]))

            # add bending moment due to element loads
            for element_load in element_loads:
                bmd[i] += element_load.get_internal_bmd(xi)

        return (stations, bmd)

    def get_sf(self, xis, analysis_case):
        """Returns the shear force within the element at *xis* for an analysis_case.

        :param xis: Positions along the element to calculate the shear force
        :type xis: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Shear forces
        :rtype: :class:`numpy.ndarray`
        """

        # ensure input is an array
        if type(xis) is float:
            xis = np.array([xis])

        # allocate shear forces
        sfs = np.zeros(len(xis))

        # get internal forces
        f = self.get_internal_actions(analysis_case=analysis_case)
        V1 = f[1]
        V2 = -f[4]

        # get list of element loads
        element_loads = self.get_element_loads(analysis_case=analysis_case)

        for (i, xi) in enumerate(xis):
            # get shape functions at xi
            (N, _) = self.get_shape_function(xi)

            # compute shear force diagram
            sf = np.dot(N, np.array([V1, V2]))

            # add shear force due to element loads
            for element_load in element_loads:
                sf += element_load.get_internal_sfd(xi)

            sfs[i] = sf

        return sfs

    def get_bm(self, xis, analysis_case):
        """Returns the bending moment within the element at *xis* for an analysis_case.

        :param xis: Positions along the element to calculate the bending moment
        :type xis: :class:`numpy.ndarray`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Bending moments
        :rtype: :class:`numpy.ndarray`
        """

        # ensure input is an array
        if type(xis) is float:
            xis = np.array([xis])

        # allocate bending moments
        bms = np.zeros(len(xis))

        # get internal forces
        f = self.get_internal_actions(analysis_case=analysis_case)
        M1 = f[2]
        M2 = -f[5]

        # get list of element loads
        element_loads = self.get_element_loads(analysis_case=analysis_case)

        for (i, xi) in enumerate(xis):
            # get shape functions at xi
            (N, _) = self.get_shape_function(xi)

            # compute bending moment diagram
            bm = np.dot(N, np.array([M1, M2]))

            # add bending moment due to element loads
            for element_load in element_loads:
                bm += element_load.get_internal_bmd(xi)

            bms[i] = bm

        return bms

    def calculate_local_displacement(self, xi, u_el):
        """Calculates the local displacement of the element at position *xi* given the displacement
        vector *u_el*.

        :param float xi: Position along the element
        :param u_el: Element displacement vector
        :type u_el: :class:`numpy.ndarray`

        :returns: Local displacement of the element *(u, v, w)*
        :rtype: tuple(float, float, float)
        """

        # element shape functions
        (N_u, N_v) = self.get_shape_function(xi)

        # compute local displacements
        u = np.dot(N_u, np.array([u_el[0, 0], u_el[1, 0]]))
        v = np.dot(N_v, np.array([u_el[0, 1], u_el[0, 2], u_el[1, 1], u_el[1, 2]]))

        return (u, v, None)

    def generate_udl(self, q):
        """Returns a EulerBernoulli2D_2N UniformDistributedLoad object for the current element.

        :param float q: Value of the uniformly distributed load

        :returns: UniformDistributedLoad object
        :rtype: :class:`~feastruct.fea.frame.EulerBernoulli2D_2N.UniformDistributedLoad`
        """

        return self.UniformDistributedLoad(self, q)

    class UniformDistributedLoad(ElementLoad):
        """Class for the application of a uniformly distributed load to a EulerBernoulli2D_2N
        element.

        :cvar element: EulerBernoulli2D_2N element to which the load is applied
        :vartype element: :class:`~feastruct.fea.frame.EulerBernoulli2D_2N`
        :cvar float q: Value of the uniformly distributed load
        """

        def __init__(self, element, q):
            """Inits the UniformDistributedLoad class.

            :param element: EulerBernoulli2D_2N element to which the load is applied
            :type element: :class:`~feastruct.fea.frame.EulerBernoulli2D_2N`
            :param float q: Value of the uniformly distributed load
            """

            super().__init__(element)
            self.q = q

        def nodal_equivalent_loads(self):
            """a"""

            # get relevant properties
            (_, _, l0, _) = self.element.get_geometric_properties()

            f_eq = np.array([
                0,
                -self.q * l0 / 2,
                -self.q * l0 * l0 / 12,
                0,
                -self.q * l0 / 2,
                self.q * l0 * l0 / 12
            ])

            return f_eq

        def apply_load(self, f_eq):
            """a"""

            # get gdofs for the element
            gdofs = self.element.get_gdof_nums()

            # calculate the nodal equivalent loads
            f_e_eq = self.nodal_equivalent_loads()

            # get relevant properties
            (_, _, _, c) = self.element.get_geometric_properties()
            cx = c[0]
            cy = c[1]

            # rotate
            f_e_eq = np.array([
                f_e_eq[0] * cx + f_e_eq[1] * cy,
                -f_e_eq[0] * cy + f_e_eq[1] * cx,
                f_e_eq[2],
                f_e_eq[3] * cx + f_e_eq[4] * cy,
                -f_e_eq[3] * cy + f_e_eq[4] * cx,
                f_e_eq[5]
            ])

            # apply fixed end forces
            f_eq[gdofs] += f_e_eq

        def get_internal_bmd(self, xi):
            """a"""

            # get relevant properties
            (_, _, l0, _) = self.element.get_geometric_properties()

            return -1 * (xi - 1) * (xi + 1) * self.q * l0 * l0 / 8

        def get_internal_sfd(self, xi):
            """a"""

            return 0

        def plot_load(self):
            """a"""

            pass
