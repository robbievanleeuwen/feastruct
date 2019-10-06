import numpy as np
from feastruct.fea.elements.frame import FrameElement


class FrameElement3D(FrameElement):
    """Class for a 3D frame element.

    Provides a number of methods that can be used for a 3D frame element.

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
        """Inits the FrameElement3D class.

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


class Bar3D_2N(FrameElement3D):
    """Two noded, three dimensional bar element that can resist an axial force only. The element
    is defined by its two end nodes and uses two linear shape functions to obtain analytical
    results.

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
        """Inits the Bar3D_2N class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # set the element freedom signature
        efs = [True, True, True, False, False, False]

        # initialise parent FrameElement3D class
        super().__init__(nodes=nodes, material=material, efs=efs, section=section)

    def get_shape_function(self, xi):
        """Returns the value of the shape functions *N1* and *N2* at *xi*.

        :param float xi: Position along the element

        :returns: Value of the shape functions *(N1, N2)* at *xi*
        :rtype: :class:`numpy.ndarray`
        """

        return np.array([0.5 - xi / 2, 0.5 + xi / 2])

    def get_stiffness_matrix(self):
        """Gets the stiffness matrix for a two noded, 3D bar element. The stiffness matrix has been
        analytically integrated so numerical integration is not necessary.

        :returns: 6 x 6 element stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, c) = self.get_geometric_properties()

        # extract relevant properties
        E = self.material.elastic_modulus
        A = self.section.area
        cx = c[0]
        cy = c[1]
        cz = c[2]

        # construct rotation matrix
        T = np.array([
            [cx, cy, cz, 0, 0, 0],
            [0, 0, 0, cx, cy, cz]
        ])

        # compute bar stiffness matrix
        k = E * A / l0 * np.array([
            [1, -1],
            [-1, 1]
        ])

        return np.matmul(np.matmul(np.transpose(T), k), T)

    def get_geometric_stiff_matrix(self, analysis_case):
        """Gets the geometric stiffness matrix for a two noded, 3D bar element. The stiffness
        matrix has been analytically integrated so numerical integration is not necessary. The
        geometric stiffness matrix requires an axial force so the analysis_case from a static
        analysis must be provided.

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
        cz = c[2]

        # get axial force
        f_int = self.get_fint(analysis_case)

        # get axial force in element (take average of nodal values)
        N = np.mean([-f_int[0], f_int[1]])

        # construct rotation matrix
        T = np.array([
            [0, 0, 0, cx, cy, cz],
            [cx, cy, cz, 0, 0, 0]
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
        cz = c[2]

        # construct rotation matrix
        T = np.array([
            [cx, cy, cz, 0, 0, 0],
            [0, 0, 0, cx, cy, cz]
        ])

        # compute element mass matrix
        m = rho * A * l0 / 6 * np.array([
            [2, 1],
            [1, 2]
        ])

        return np.matmul(np.matmul(np.transpose(T), m), T)

    def get_internal_actions(self, analysis_case):
        """Returns the internal actions for a two noded 3D bar element.

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
        cz = c[2]

        f = np.array([
            f_int[0] * cx + f_int[1] * cy + f_int[2] * cz,
            f_int[3] * cx + f_int[4] * cy + f_int[5] * cz
        ])

        return f

    # def calculate_local_displacement(self, xi, u_el, analysis_case):
    #     """Calculates the local displacement of the element at position *xi* given the displacement
    #     vector *u_el* for a Bar3D_2N element.
    #
    #     :param float xi: Position along the element
    #     :param u_el: Element displacement vector
    #     :type u_el: :class:`numpy.ndarray`
    #     :param analysis_case: Analysis case
    #     :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    #
    #     :returns: Local displacement of the element *(u, v, w, rx, ry, rz)*
    #     :rtype: tuple(float)
    #     """
    #
    #     # element shape function
    #     N = self.get_shape_function(xi)
    #
    #     # compute local displacements
    #     u = np.dot(N, np.array([u_el[0, 0], u_el[1, 0]]))
    #     v = np.dot(N, np.array([u_el[0, 1], u_el[1, 1]]))
    #     w = np.dot(N, np.array([u_el[0, 2], u_el[1, 2]]))
    #
    #     return (u, v, w, None, None, None)

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


class EulerBernoulli3D_2N(FrameElement3D):
    """Two noded, three dimensional frame element based on the Euler-Bernoulli beam formulation for
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
        """Inits the EulerBernoulli3D_2N class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # set the element freedom signature
        efs = [True, True, True, True, True, True]

        # initialise parent FrameElement3D class
        super().__init__(nodes=nodes, material=material, efs=efs, section=section)

    def get_shape_function(self, xi):
        """Returns the value of the shape functions *Nu1*, *Nu2*, *Nv1* to *Nv4* at *xi*.

        :param float xi: Position along the element

        :returns: Value of the shape functions *((Nu1, Nu2), (Nv1 to Nv4))* at *xi*
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

    def get_stiffness_matrix(self, linear=True):
        """Gets the stiffness matrix for a two noded 3D Euler-Bernoulli frame element. The
        stiffness matrix has been analytically integrated so numerical integration is not
        necessary.

        :returns: 12 x 12 element stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, _) = self.get_geometric_properties()

        # extract relevant properties
        E = self.material.elastic_modulus
        G = self.material.shear_modulus
        A = self.section.area
        ixx = self.section.ixx
        iyy = self.section.iyy
        j = self.section.j
        l02 = l0 * l0
        l03 = l0 * l0 * l0

        # compute stiffness matrix
        k_el = np.array([
            [E*A/l0, 0, 0, 0, 0, 0, -E*A/l0, 0, 0, 0, 0, 0],
            [0, 12*E*ixx/l03, 0, 0, 0, 6*E*ixx/l02, 0, -12*E*ixx/l03, 0, 0, 0, 6*E*ixx/l02],
            [0, 0, 12*E*iyy/l03, 0, -6*E*iyy/l02, 0, 0, 0, -12*E*iyy/l03, 0, -6*E*iyy/l02, 0],
            [0, 0, 0, G*j/l0, 0, 0, 0, 0, 0, -G*j/l0, 0, 0],
            [0, 0, -6*E*iyy/l02, 0, 4*E*iyy/l0, 0, 0, 0, 6*E*iyy/l02, 0, 2*E*iyy/l0, 0],
            [0, 6*E*ixx/l02, 0, 0, 0, 4*E*ixx/l0, 0, -6*E*ixx/l02, 0, 0, 0, 2*E*ixx/l0],
            [-E*A/l0, 0, 0, 0, 0, 0, E*A/l0, 0, 0, 0, 0, 0],
            [0, -12*E*ixx/l03, 0, 0, 0, -6*E*ixx/l02, 0, 12*E*ixx/l03, 0, 0, 0, -6*E*ixx/l02],
            [0, 0, -12*E*iyy/l03, 0, 6*E*iyy/l02, 0, 0, 0, 12*E*iyy/l03, 0, 6*E*iyy/l02, 0],
            [0, 0, 0, -G*j/l0, 0, 0, 0, 0, 0, G*j/l0, 0, 0],
            [0, 0, -6*E*iyy/l02, 0, 2*E*iyy/l0, 0, 0, 0, 6*E*iyy/l02, 0, 4*E*iyy/l0, 0],
            [0, 6*E*ixx/l02, 0, 0, 0, 2*E*ixx/l0, 0, -6*E*ixx/l02, 0, 0, 0, 4*E*ixx/l0]
        ])

        # construct rotation matrix
        T = self.get_transformation_matrix()

        return np.matmul(np.matmul(np.transpose(T), k_el), T)

    def get_geometric_stiff_matrix(self, analysis_case):
        """Gets the geometric stiffness matrix for a two noded 3D Euler-Bernoulli frame element.
        The stiffness matrix has been analytically integrated so numerical integration is not
        necessary. The geometric stiffness matrix requires an axial force so the analysis_case from
        a static analysis must be provided.

        :param analysis_case: Analysis case from which to extract the axial force
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :returns: 12 x 12 element geometric stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        # compute geometric parameters
        (_, _, l0, _) = self.get_geometric_properties()

        # extract relevant properties
        A = self.section.area
        j = self.section.j

        # get axial force
        f_int = self.get_fint(analysis_case)

        # get axial force in element (take average of nodal values)
        N = np.mean([-f_int[0], f_int[6]])

        # form geometric stiffness matrix
        k_el_g = np.array([
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1.2, 0, 0, 0, 0.1*l0, 0, -1.2, 0, 0, 0, 0.1*l0],
            [0, 0, 1.2, 0, -0.1*l0, 0, 0, 0, -1.2, 0, -0.1*l0, 0],
            [0, 0, 0, j/A, 0, 0, 0, 0, 0, -j/A, 0, 0],
            [0, 0, -0.1*l0, 0, 2*l0*l0/15, 0, 0, 0, 0.1*l0, 0, -l0*l0/30, 0],
            [0, 0.1*l0, 0, 0, 0, 2*l0*l0/15, 0, -0.1*l0, 0, 0, 0, -l0*l0/30],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, -1.2, 0, 0, 0, -0.1*l0, 0, 1.2, 0, 0, 0, -0.1*l0],
            [0, 0, -1.2, 0, 0.1*l0, 0, 0, 0, 1.2, 0, 0.1*l0, 0],
            [0, 0, 0, -j/A, 0, 0, 0, 0, 0, j/A, 0, 0],
            [0, 0, -0.1*l0, 0, -l0*l0/30, 0, 0, 0, 0.1*l0, 0, 2*l0*l0/15, 0],
            [0, 0.1*l0, 0, 0, 0, -l0*l0/30, 0, -0.1*l0, 0, 0, 0, 2*l0*l0/15]
        ])
        k_el_g *= N / l0

        # construct rotation matrix
        T = self.get_transformation_matrix()

        return np.matmul(np.matmul(np.transpose(T), k_el_g), T)

    # def get_mass_matrix(self):
    #     """Gets the mass matrix for a for a two noded 2D Euler-Bernoulli frame element. The mass
    #     matrix has been analytically integrated so numerical integration is not necessary.
    #
    #     :returns: 6 x 6 element mass matrix
    #     :rtype: :class:`numpy.ndarray`
    #     """
    #
    #     # compute geometric parameters
    #     (_, _, l0, c) = self.get_geometric_properties()
    #
    #     # extract relevant properties
    #     rho = self.material.rho
    #     A = self.section.area
    #     cx = c[0]
    #     cy = c[1]
    #
    #     # compute element mass matrix
    #     m_el = np.array([
    #         [140, 0, 0, 70, 0, 0],
    #         [0, 156, 22*l0, 0, 54, -13*l0],
    #         [0, 22*l0, 4*l0*l0, 0, 13*l0, -3*l0*l0],
    #         [70, 0, 0, 140, 0, 0],
    #         [0, 54, 13*l0, 0, 156, -22*l0],
    #         [0, -13*l0, -3*l0*l0, 0, -22*l0, 4*l0*l0]
    #     ])
    #     m_el *= rho * A * l0 / 420
    #
    #     # construct rotation matrix
    #     T = np.array([
    #         [cx, cy, 0, 0, 0, 0],
    #         [-cy, cx, 0, 0, 0, 0],
    #         [0, 0, 1, 0, 0, 0],
    #         [0, 0, 0, cx, cy, 0],
    #         [0, 0, 0, -cy, cx, 0],
    #         [0, 0, 0, 0, 0, 1]
    #     ])
    #
    #     return np.matmul(np.matmul(np.transpose(T), m_el), T)
    #
    # def get_internal_actions(self, analysis_case):
    #     """Returns the internal actions for a two noded 3D Euler-Bernoulli frame element.
    #
    #     :param analysis_case: Analysis case
    #     :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    #
    #     :returns: An array containing the internal actions for the element
    #         *(N1, V1, M1, N2, V2, M2)*
    #     :rtype: :class:`numpy.ndarray`
    #     """
    #
    #     (_, _, _, c) = self.get_geometric_properties()
    #     f_int = self.get_fint(analysis_case=analysis_case)
    #
    #     cx = c[0]
    #     cy = c[1]
    #
    #     f = np.array([
    #         f_int[0] * cx + f_int[1] * cy,
    #         -f_int[0] * cy + f_int[1] * cx,
    #         f_int[2],
    #         f_int[3] * cx + f_int[4] * cy,
    #         -f_int[3] * cy + f_int[4] * cx,
    #         f_int[5]
    #     ])
    #
    #     return f
    #
    # def calculate_local_displacement(self, xi, u_el, analysis_case):
    #     """Calculates the local displacement of the element at position *xi* given the displacement
    #     vector *u_el*.
    #
    #     :param float xi: Position along the element
    #     :param u_el: Element displacement vector
    #     :type u_el: :class:`numpy.ndarray`
    #     :param analysis_case: Analysis case
    #     :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    #
    #     :returns: Local displacement of the element *(u, v)*
    #     :rtype: tuple(float, float)
    #     """
    #
    #     # element shape functions
    #     (N_u, N_v) = self.get_shape_function(xi)
    #
    #     # compute local displacements
    #     u = np.dot(N_u, np.array([u_el[0, 0], u_el[1, 0]]))
    #     v = np.dot(N_v, np.array([u_el[0, 1], u_el[0, 2], u_el[1, 1], u_el[1, 2]]))
    #
    #     return (u, v)

    def get_transformation_matrix(self):
        """Returns the transformation matrix for an EulerBernoulli3D_2N element.

        :returns: Element transformation matrix
        :rtype: :class:`numpy.ndarray`
        """

        (_, _, _, c) = self.get_geometric_properties()

        # construct rotation matrix
        cxx = c[0]
        cyx = c[1]
        czx = c[2]
        d = np.sqrt(cxx * cxx + cyx * cyx)
        cxy = -cyx / d
        cyy = cxx / d
        czy = 0
        cxz = -cxx * czx / d
        cyz = -cyx * czx / d
        czz = d

        gamma = np.array([
            [cxx, cyx, czx],
            [cxy, cyy, czy],
            [cxz, cyz, czz]
        ])

        T = np.zeros((12, 12))
        T[0:3, 0:3] = gamma
        T[3:6, 3:6] = gamma
        T[6:9, 6:9] = gamma
        T[9:12, 9:12] = gamma

        return T

    #
    # def generate_udl(self, q):
    #     """Returns a EulerBernoulli2D_2N UniformDistributedLoad object for the current element.
    #
    #     :param float q: Value of the uniformly distributed load
    #
    #     :returns: UniformDistributedLoad object
    #     :rtype: :class:`~feastruct.fea.frame.EulerBernoulli2D_2N.UniformDistributedLoad`
    #     """
    #
    #     return self.UniformDistributedLoad(self, q)
    #
    # class UniformDistributedLoad(ElementLoad):
    #     """Class for the application of a uniformly distributed load to a EulerBernoulli2D_2N
    #     element.
    #
    #     :cvar element: EulerBernoulli2D_2N element to which the load is applied
    #     :vartype element: :class:`~feastruct.fea.frame.EulerBernoulli2D_2N`
    #     :cvar float q: Value of the uniformly distributed load
    #     """
    #
    #     def __init__(self, element, q):
    #         """Inits the UniformDistributedLoad class.
    #
    #         :param element: EulerBernoulli2D_2N element to which the load is applied
    #         :type element: :class:`~feastruct.fea.frame.EulerBernoulli2D_2N`
    #         :param float q: Value of the uniformly distributed load
    #         """
    #
    #         super().__init__(element)
    #         self.q = q
    #
    #     def nodal_equivalent_loads(self):
    #         """a"""
    #
    #         # get relevant properties
    #         (_, _, l0, _) = self.element.get_geometric_properties()
    #
    #         f_eq = np.array([
    #             0,
    #             -self.q * l0 / 2,
    #             -self.q * l0 * l0 / 12,
    #             0,
    #             -self.q * l0 / 2,
    #             self.q * l0 * l0 / 12
    #         ])
    #
    #         return f_eq
    #
    #     def apply_load(self, f_eq):
    #         """a"""
    #
    #         # get gdofs for the element
    #         gdofs = self.element.get_gdof_nums()
    #
    #         # calculate the nodal equivalent loads
    #         f_e_eq = self.nodal_equivalent_loads()
    #
    #         # get relevant properties
    #         (_, _, _, c) = self.element.get_geometric_properties()
    #         cx = c[0]
    #         cy = c[1]
    #
    #         # rotate
    #         f_e_eq = np.array([
    #             f_e_eq[0] * cx + f_e_eq[1] * cy,
    #             -f_e_eq[0] * cy + f_e_eq[1] * cx,
    #             f_e_eq[2],
    #             f_e_eq[3] * cx + f_e_eq[4] * cy,
    #             -f_e_eq[3] * cy + f_e_eq[4] * cx,
    #             f_e_eq[5]
    #         ])
    #
    #         # apply fixed end forces
    #         f_eq[gdofs] += f_e_eq
    #
    #     def get_internal_bmd(self, xi):
    #         """a"""
    #
    #         # get relevant properties
    #         (_, _, l0, _) = self.element.get_geometric_properties()
    #
    #         return -1 * (xi - 1) * (xi + 1) * self.q * l0 * l0 / 8
    #
    #     def get_internal_sfd(self, xi):
    #         """a"""
    #
    #         return 0
    #
    #     def plot_load(self):
    #         """a"""
    #
    #         pass
