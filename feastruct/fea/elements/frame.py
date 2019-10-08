import numpy as np
from scipy import optimize
from feastruct.fea.fea import FiniteElement


class FrameElement(FiniteElement):
    """Parent class for a frame element.

    Establishes a template for a frame element and provides a number of generic methods that can be
    used for any frame element.

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

    def __init__(self, nodes, material, efs, section):
        """Inits the FrameElement class.

        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param efs: Element freedom signature
        :type efs: list[bool]
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        """

        # initialise parent FiniteElement class
        super().__init__(nodes=nodes, material=material, efs=efs)

        self.section = section

    def get_geometric_properties(self):
        """Calculates geometric properties related to a frame element. Returns the following:

        * *node_coords*: An *(n x 3)* array of node coordinates, where n is the number of nodes for
          the given finite element
        * *dx*: A *(1 x 3)* array consisting of the x, y and z distances between the nodes
        * *l0*: The original length of the frame element
        * *c*: A *(1 x 3)* array consisting of the cosines with respect to the x, y and z axes

        :returns: *(node_coords, dx, l0, c)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`, float,
            :class:`numpy.ndarray`)
        """

        node_coords = self.get_node_coords()
        dx = node_coords[1] - node_coords[0]
        l0 = np.linalg.norm(dx)
        c = dx / l0

        return (node_coords, dx, l0, c)

    def get_sampling_points(self, n, analysis_case, bm=False, defl=False):
        """Returns a list of sampling points along a 2D frame element given a minimum *n* points
        and an analysis case. The sampling points vary from 0 to 1.

        Adds a sampling point at the following locations:

        * Concentrated element load
        * Point of zero shear force (if bm=True)
        * Point of zero rotation (if defl=True)

        :param int n: Minimum number of sampling points
        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param bool bm: Whether or not the sampling points are for the bending moment (i.e. include
            sfd roots)
        :param bool defl: Whether or not the sampling points for transverse deflection (i.e.
            include rotation roots). N.B. the FrameElement must have a `calculate_rotation` method.

        :returns: List of sampling points
        :rtype: :class:`numpy.ndarray`
        """

        # generate initial list of stations
        stations = np.linspace(0, 1, n)

        # find any points of zero shear force
        if bm:
            # get sfd
            (xis, sfd) = self.get_sfd(n=n, analysis_case=analysis_case)

            # loop through shear force diagram
            for i in range(len(xis) - 1):
                # if there is a root between i and i + 1 (different signs)
                if (sfd[i] > 0 and sfd[i+1] < 0) or (sfd[i] < 0 and sfd[i+1] > 0):
                    # determine root using brentq method
                    def sf(x): return self.get_sf(x, analysis_case)

                    # search for root between two points
                    xi = optimize.brentq(sf, xis[i], xis[i+1])

                    # check that the station doesn't already exist
                    for diff in (stations - xi):
                        # if the station already exists
                        if abs(diff) < 1e-5:
                            break
                    else:
                        # add the station if it doesn't already exist
                        stations = np.append(stations, xi)

        # if defl - find any points of zero rotation
        if defl:
            # get the nodal displacements
            u_el = self.get_nodal_displacements(analysis_case)

            # get initial rotation
            phi0 = u_el[0, 2]

            # get rotations
            rots = self.calculate_rotation(xis=stations, phi0=phi0, analysis_case=analysis_case)

            # loop through rotations
            for i in range(len(stations) - 1):
                # if there is a root between i and i + 1 (different signs)
                if (rots[i] > 0 and rots[i+1] < 0) or (rots[i] < 0 and rots[i+1] > 0):
                    # determine root using brentq method
                    def rot(x): return self.calculate_rotation(x, phi0, analysis_case)

                    # search for root between two points
                    xi = optimize.brentq(rot, stations[i], stations[i+1])

                    # check that the station doesn't already exist
                    for diff in (stations - xi):
                        # if the station already exists
                        if abs(diff) < 1e-5:
                            break
                    else:
                        # add the station if it doesn't already exist
                        stations = np.append(stations, xi)

        # re-sort stations list
        return np.sort(stations)

    def get_element_loads(self, analysis_case):
        """Returns a list of element loads on a FrameElement for an analyis case.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        element_loads = []  # list of applied element loads

        for element_load in analysis_case.load_case.element_items:
            # if the current element has an applied element load
            if element_load.element is self:
                # add nodal equivalent loads to f_int
                element_loads.append(element_load)

        return element_loads

    def get_displacements(self, n, analysis_case):
        """Placeholder for the get_displacements method.

        Returns a list of the local displacements, *(u, v, w, ru, rv, rw)*, along the element for
        the analysis case and a minimum of *n* subdivisions. A list of the stations, *xi*, is also
        included. Station locations, *xis*, vary from 0 to 1.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int n: Minimum number of sampling points

        :returns: 2D numpy array containing stations and local displacements. Each station contains
            an array of the following format: *[xi, u, v, w, rx, ry, rz]*
        :rtype: :class:`numpy.ndarray`
        """

        pass

    def get_transformation_matrix(self):
        """Placeholder for the get_transformation_matrix method.

        Returns the transformation matrix for a FrameElement.

        :returns: Element transformation matrix
        :rtype: :class:`numpy.ndarray`
        """

        pass

    def get_afd(self, n, analysis_case):
        """Placeholder for the get_afd method.

        Returns the axial force diagram within the element for a minimum of *n* stations for an
        analysis_case. Station locations, *xis*, vary from 0 to 1.

        :param int n: Minimum number of stations to sample the axial force diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and axial force diagram, *afd* - *(xis, afd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        pass

    def get_sfd(self, n, analysis_case):
        """Placeholder for the get_sfd method.

        Returns the shear force diagram within the element for a minimum of *n* stations for an
        analysis_case. Station locations, *xis*, vary from 0 to 1.

        :param int n: Minimum number of stations to sample the shear force diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and shear force diagram, *sfd* - *(xis, sfd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        pass

    def get_bmd(self, n, analysis_case):
        """Placeholder for the get_bmd method.

        Returns the bending moment diagram within the element for a minimum of *n* stations for
        an analysis_case. Station locations, *xis*, vary from 0 to 1.

        :param int n: Minimum number of stations to sample the bending moment diagram
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Station locations, *xis*, and bending moment diagram, *bmd* - *(xis, bmd)*
        :rtype: tuple(:class:`numpy.ndarray`, :class:`numpy.ndarray`)
        """

        pass

    def calculate_local_displacement(self, xi, u_el):
        """Placeholder for the calculate_local_displacement method.

        Calculates the local displacement of the element at position *xi* given the displacement
        vector *u_el*.

        :param float xi: Position along the element *(0 < xi < 1)*
        :param u_el: Element displacement vector
        :type u_el: :class:`numpy.ndarray`

        :returns: Local displacement of the element *(u, v, w)*
        :rtype: tuple(float, float, float)
        """

        pass

    def map_to_station(self, eta):
        """Maps the isometric parameter *-1 < eta < 1* to a station value *0 < xi < 1*.

        :param float xi: Isoparametric coordinate (*-1 < eta < 1*)

        :returns: Station location (*0 < x < 1*)
        :rtype: float
        """

        return 0.5 * (eta + 1)

    def map_to_isoparam(self, xi):
        """Maps a station value *0 < xi < 1* to the isometric parameter *-1 < eta < 1*.

        :param float xi: Station location (*0 < x < 1*)

        :returns: Isoparametric coordinate (*-1 < eta < 1*)
        :rtype: float
        """

        return 2 * xi - 1
