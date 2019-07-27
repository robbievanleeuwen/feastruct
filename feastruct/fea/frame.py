import numpy as np
from matplotlib.patches import Polygon
from feastruct.fea.fea import FiniteElementAnalysis, FiniteElement


class FrameAnalysis(FiniteElementAnalysis):
    """Parent class for a frame analysis.

    Includes a method for analysis initiation and a method for adding a frame element to the
    analysis.

    :cvar nodes: Nodes used in the finite element analysis
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis
    :vartype elements: list[:class:`~feastruct.fea.frame.FrameElement`]
    :cvar int dims: Number of dimensions used for the current analysis type
    :cvar dofs: List of the degrees of freedom used in the current analysis type
    :vartype dofs: list[int]
    :cvar post: Post-processor object
    :vartype post: :class:`feastruct.post.post.PostProcessor`
    """

    def __init__(self, nodes, elements, dims, dofs):
        """Inits the FrameAnalysis class.

        :param nodes: List of nodes with which to initialise the class
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class
        :type elements: list[:class:`~feastruct.fea.frame.FrameElement`]
        :param int dims: Number of dimensions used for the current analysis type
        :param dofs: List of the degrees of freedom used in the current analysis type
        :type dofs: list[int]
        """

        # initialise parent FiniteElementAnalysis class
        super().__init__(nodes=nodes, elements=elements, dims=dims, dofs=dofs)

    def create_element(self, el_type, nodes, material, section):
        """Creates and returns a frame element and adds it to the
        :class:`~feastruct.fea.frame.Frame` object. Refer to 'xxx' for the possible element types.

        :param string el_type: String characterising the type of frame element
        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        :return: FrameElement object
        :rtype: :class:`~feastruct.fea.frame.FrameElement`
        """

        if el_type == 'Bar2-2D':
            element = Bar2D_2N(nodes=nodes, material=material, section=section)
        elif el_type == 'EB2-2D':
            element = EulerBernoulli2D_2N(nodes=nodes, material=material, section=section)

        return(FiniteElementAnalysis.create_element(self, element=element))


class FrameAnalysis2D(FrameAnalysis):
    """Parent class for a 2D frame analysis.

    Includes a method for analysis initiation and a method for adding 2D frame elements to the
    analysis.

    :cvar nodes: Nodes used in the finite element analysis
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis
    :vartype elements: list[:class:`~feastruct.fea.frame.FrameElement`]
    :cvar int dims: Number of dimensions used for the current analysis type
    :cvar dofs: List of the degrees of freedom used in the current analysis type
    :vartype dofs: list[int]
    :cvar post: Post-processor object
    :vartype post: :class:`feastruct.post.post.PostProcessor`
    """

    def __init__(self, nodes=None, elements=None):
        """Inits the FrameAnalysis2D class.

        :param nodes: List of nodes with which to initialise the class
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class
        :type elements: list[:class:`~feastruct.fea.frame.FrameElement`]
        """

        # initialise parent FrameAnalysis class
        super().__init__(nodes=nodes, elements=elements, dims=2, dofs=[0, 1, 5])

    def create_element(self, el_type, nodes, material, section):
        """Creates and returns a 2D frame element and adds it to the
        :class:`~feastruct.fea.frame.Frame2D` object. Refer to 'xxx' for the possible 2D element
        types.

        :param string el_type: String characterising the type of frame element
        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        :return: FrameElement object
        :rtype: :class:`~feastruct.fea.frame.FrameElement`
        """

        # TODO: check to ensure element is 2D

        return(FrameAnalysis.create_element(
            self, el_type=el_type, nodes=nodes, material=material, section=section
        ))


class FrameAnalysis3D(FrameAnalysis):
    pass


class FrameElement(FiniteElement):
    """Parent class for a frame element.

    Establishes a template for a frame element and provides a number of generic methods that can be
    used for any frame element.

    :cvar nodes: List of node objects defining the element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
    :cvar f_int: List of internal force vector results stored for each analysis case
    :vartype f_int: list[:class:`~feastruct.fea.fea.ForceVector`]
    :cvar section: Section object for the element
    :vartype section: :class:`~feastruct.pre.section.Section`
    """

    def __init__(self, nodes, material, section):
        """Inits the FrameElement class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # initialise parent FiniteElement class
        super().__init__(nodes=nodes, material=material)

        self.section = section

    def get_geometric_properties(self):
        """Calculates geometric properties related to a frame element. Returns
        the following:

        * *node_coords*: An *(n x 3)* array of node coordinates, where n is the number of nodes for
          the given finite element
        * *dx*: A *(1 x 3)* array consisting of the x, y and z distances between the nodes
        * *l0*: The original length of the frame element
        * *c*: A *(1 x 3)* array consisting of the cosines with respect to the x, y and z axes

        :return: *(node_coords, dx, l0, c)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, float,
            :class:`numpy.ndarray`)
        """

        node_coords = self.get_node_coords()
        dx = node_coords[1] - node_coords[0]
        l0 = np.linalg.norm(dx)
        c = dx / l0

        return (node_coords, dx, l0, c)


class FrameElement2D(FrameElement):
    """Class for a 2D frame element.

    Provides a number of methods that can be used for a 2D frame element.

    :cvar nodes: List of node objects defining the element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
    :cvar f_int: List of internal force vector results stored for each analysis
        case
    :vartype f_int: list[:class:`~feastruct.fea.fea.ForceVector`]
    :cvar section: Section object for the element
    :vartype section: :class:`~feastruct.pre.section.Section`
    """

    def __init(self, nodes, material, section):
        """Inits the FrameElement2D class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # initialise parent FrameElement class
        super().__init__(nodes=nodes, material=material, section=section)

    def get_nodal_displacements(self, analysis_case):
        """Returns an array of the nodal displacements for each degree of freedom in the 2D frame
        element *([0, 1, 5])* for the analysis_case.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :return: An *(n_nodes x 3)* array of degree of freedom displacements
        :rtype: :class:`numpy.ndarray`
        """

        # get nodal displacements for x & y translations and rotation about z
        return(super().get_nodal_displacements(dof_nums=[0, 1, 5], analysis_case=analysis_case))

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

    def plot_deformed_element(self, ax, analysis_case, n, def_scale):
        """Plots a 2D frame element in its deformed configuration for the displacement vector
        defined by the analysis_case. The deformation is based on the element shape functions.

        :param ax: Axis object on which to draw the element
        :type ax: :class:`matplotlib.axes.Axes`
        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int n: Number of linear subdivisions used to plot the element
        :param float def_scale: Deformation scale
        """

        # TODO: IMPLEMENT!
        # get nodal displacements
        u_el = self.get_nodal_displacements(analysis_case)

        # compute frame geometric parameters
        (node_coords, _, l0, c) = self.get_geometric_properties()

        # rotate nodal displacements to local axis
        T = np.array([[c[0], c[1], 0], [-c[1], c[0], 0], [0, 0, 1]])
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
            N_v = np.array([
                0.25 * (1 - xi) * (1 - xi) * (2 + xi),
                0.125 * l0 * (1 - xi) * (1 - xi) * (1 + xi),
                0.25 * (1 + xi) * (1 + xi) * (2 - xi),
                -0.125 * l0 * (1 + xi) * (1 + xi) * (1 - xi)
            ])

            # compute local displacements
            u = np.dot(N_u, np.array([u_el[0, 0], u_el[1, 0]]))
            v = np.dot(N_v, np.array([u_el[0, 1], u_el[0, 2],
                                      u_el[1, 1], u_el[1, 2]]))

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

    def plot_axial_force(self, ax, analysis_case, scalef):
        """Plots the axial force diagram from a static analysis defined by case_id. N.B. this
        method is adopted from the MATLAB code by F.P. van der Meer: plotNLine.m.

        :param ax: Axis object on which to draw the element
        :type ax: : class: `matplotlib.axes.Axes`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param float scalef: Factor by which to scale the axial force diagram
        """

        # get geometric properties
        (node_coords, dx, l0, _) = self.get_geometric_properties()

        # get internal force vector
        f_int = self.get_internal_actions(analysis_case=analysis_case)
        n1 = -f_int[0]  # axial force at node 1 (tension positive)
        n2 = f_int[3]  # axial force at node 2 (tension positive)

        # location of node 1 and node 2
        p1 = node_coords[0, 0:2]
        p2 = node_coords[1, 0:2]

        # location of the axial force diagram end points
        v = np.matmul(np.array([[0, -1], [1, 0]]), dx[0:2]) / l0  # direction vector
        p3 = p2 + v * scalef * n2
        p4 = p1 + v * scalef * n1

        # plot axial force line and patch
        ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0.7, 0, 0))
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0.7, 0, 0))
        ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0.7, 0, 0))
        ax.add_patch(
            Polygon(np.array([p1, p2, p3, p4]), facecolor=(1, 0, 0), linestyle='None', alpha=0.3))

        # plot text value of axial force
        mid1 = (p1 + p4) / 2
        mid2 = (p2 + p3) / 2
        ax.text(mid1[0], mid1[1], "{:5.3g}".format(n1), size=8, verticalalignment='bottom')
        ax.text(mid2[0], mid2[1], "{:5.3g}".format(n2), size=8, verticalalignment='bottom')

    def plot_shear_force(self, ax, analysis_case, scalef):
        """Plots the axial force diagram from a static analysis defined by case_id. N.B. this
        method is adopted from the MATLAB code by F.P. van der Meer: plotVLine.m.

        :param ax: Axis object on which to draw the element
        :type ax: : class: `matplotlib.axes.Axes`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param float scalef: Factor by which to scale the shear force diagram
        """

        # get geometric properties
        (node_coords, dx, l0, _) = self.get_geometric_properties()

        # get internal force vector
        f_int = self.get_internal_actions(analysis_case=analysis_case)
        v1 = f_int[1]  # shear force at node 1 (cw positive)
        v2 = -f_int[4]  # shear force at node 2 (cw positive)

        # location of node 1 and node 2
        p1 = node_coords[0, 0:2]
        p2 = node_coords[1, 0:2]

        # location of the shear force diagram end points
        v = np.matmul(np.array([[0, -1], [1, 0]]), dx[0:2]) / l0  # direction vector
        p3 = p2 + v * scalef * v2
        p4 = p1 + v * scalef * v1

        # plot shear force line and patch
        ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0, 0.3, 0))
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0, 0.3, 0))
        ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0, 0.3, 0))
        ax.add_patch(Polygon(
            np.array([p1, p2, p3, p4]), facecolor=(0, 0.5, 0), linestyle='None', alpha=0.3))

        # plot text value of shear force
        mid1 = (p1 + p4) / 2
        mid2 = (p2 + p3) / 2
        ax.text(mid1[0], mid1[1], "{:5.3g}".format(v1), size=8, verticalalignment='bottom')
        ax.text(mid2[0], mid2[1], "{:5.3g}".format(v2), size=8, verticalalignment='bottom')

    def plot_bending_moment(self, ax, analysis_case, scalef):
        """Plots the axial force diagram from a static analysis defined by case_id. N.B. this
        method is adopted from the MATLAB code by F.P. van der Meer: plotMLine.m.

        :param ax: Axis object on which to draw the element
        :type ax: : class: `matplotlib.axes.Axes`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param float scalef: Factor by which to scale the bending moment diagram
        """

        # get geometric properties
        (node_coords, dx, l0, _) = self.get_geometric_properties()

        # get internal force vector
        f_int = self.get_internal_actions(analysis_case=analysis_case)
        m1 = f_int[2]  # bending moment at node 1 (tension positive)
        m2 = -f_int[5]  # bending moment at node 2 (tension positive)

        # location of node 1 and node 2
        p1 = node_coords[0, 0:2]
        p2 = node_coords[1, 0:2]

        # location of the bending moment diagram end points
        v = np.matmul(np.array([[0, -1], [1, 0]]), dx[0:2]) / l0  # direction vector
        p3 = p2 + v * scalef * m2
        p4 = p1 + v * scalef * m1

        # plot bending moment line and patch
        ax.plot([p1[0], p4[0]], [p1[1], p4[1]], linewidth=1, color=(0, 0, 0.7))
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], linewidth=1, color=(0, 0, 0.7))
        ax.plot([p3[0], p2[0]], [p3[1], p2[1]], linewidth=1, color=(0, 0, 0.7))
        ax.add_patch(Polygon(
            np.array([p1, p2, p3, p4]), facecolor=(0.2, 0.4, 0.8), linestyle='None', alpha=0.3))

        # plot text value of bending moment
        mid1 = (p1 + p4) / 2
        mid2 = (p2 + p3) / 2
        ax.text(mid1[0], mid1[1], "{:5.3g}".format(m1), size=8, verticalalignment='bottom')
        ax.text(mid2[0], mid2[1], "{:5.3g}".format(m2), size=8, verticalalignment='bottom')


class FrameElement3D(FrameElement):
    pass


class Bar2D_2N(FrameElement2D):
    """Two noded, two dimensional bar element that can resist an axial force only. The element is
    defined by its two end nodes and uses two linear shape functions to obtain analytical results.

    :cvar nodes: List of node objects defining the element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
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

        # initialise parent FrameElement2D class
        super().__init__(nodes=nodes, material=material, section=section)

    def get_stiffness_matrix(self, linear=True):
        """Gets the stiffness matrix for a two noded, 2D bar element. The stiffness matrix has been
        analytically integrated so numerical integration is not necessary.

        :param bool linear: Whether a linear or non-linear stiffness matrix is required
        :return: 4 x 4 element stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        if linear:
            # compute geometric parameters
            (_, _, l0, c) = self.get_geometric_properties()

            # extract relevant properties
            E = self.material.elastic_modulus
            A = self.section.area
            s = c[1]
            c = c[0]

            # compute bar stiffness matrix
            return E * A / l0 * np.array([
                [c*c, c*s, -c*c, -c*s],
                [c*s, s*s, -c*s, -s*s],
                [-c*c, -c*s, c*c, c*s],
                [-c*s, -s*s, c*s, s*s]
            ])
        else:
            # TODO: implement non-linear stiffness matrix
            pass

    def get_geometric_stiff_matrix(self, analysis_case):
        pass

    def get_mass_matrix(self):
        pass


class Bar3D_2N(FrameElement3D):
    pass


class EulerBernoulli2D_2N(FrameElement2D):
    """Two noded, two dimensional frame element based on the Euler-Bernoulli beam formulation for
    relatively thin beams. The element is defined by its two end nodes and uses four cubic
    polynomial shape functions to obtain analytical results.

    :cvar nodes: List of node objects defining the element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
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

        # initialise parent FrameElement2D class
        super().__init__(nodes=nodes, material=material, section=section)

    def get_stiffness_matrix(self, linear=True):
        """Gets the stiffness matrix for a two noded 2D Euler-Bernoulli frame element. The
        stiffness matrix has been analytically integrated so numerical integration is not
        necessary.

        :param bool linear: Whether a linear or non-linear stiffness matrix is required
        :return: 6 x 6 element stiffness matrix
        :rtype: :class: `numpy.ndarray`
        """

        if linear:
            # compute geometric parameters
            (_, _, l0, c) = self.get_geometric_properties()

            # extract relevant properties
            E = self.material.elastic_modulus
            A = self.section.area
            ixx = self.section.ixx
            s = c[1]
            c = c[0]

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

        else:
            pass
            # TODO: implement non-linear stiffness matrix

        # construct rotation matrix
        T = np.array([
            [c, s, 0, 0, 0, 0],
            [-s, c, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, s, 0],
            [0, 0, 0, -s, c, 0],
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
        :return: 6 x 6 element geometric stiffness matrix
        :rtype: :class: `numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        s = c[1]
        c = c[0]

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
            [c, s, 0, 0, 0, 0],
            [-s, c, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, s, 0],
            [0, 0, 0, -s, c, 0],
            [0, 0, 0, 0, 0, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), k_el_g), T)

    def get_mass_matrix(self):
        """Gets the mass matrix for a for a two noded 2D Euler-Bernoulli frame element. The mass
        matrix has been analytically integrated so numerical integration is not necessary.

        :return: 6 x 6 element mass matrix
        :rtype: : class: `numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        s = c[1]
        c = c[0]
        rho = self.material.rho
        A = self.section.area

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
            [c, s, 0, 0, 0, 0],
            [-s, c, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, c, s, 0],
            [0, 0, 0, -s, c, 0],
            [0, 0, 0, 0, 0, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), m_el), T)

    def get_internal_actions(self, analysis_case):
        """Returns the internal actions for a frame element.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: An array containing the internal actions for the element
            *(N1, V1, M1, N2, V2, M2)*
        :rtype: :class:`numpy.ndarray`
        """

        (_, _, _, c) = self.get_geometric_properties()
        f_int = self.get_fint(analysis_case=analysis_case)

        s = c[1]
        c = c[0]

        f = np.array([
            f_int[0] * c + f_int[1] * s,
            -f_int[0] * s + f_int[1] * c,
            f_int[2],
            f_int[3] * c + f_int[4] * s,
            -f_int[3] * s + f_int[4] * c,
            f_int[5]
        ])

        return f


class EulerBernoulli3D_2N(FrameElement3D):
    pass
