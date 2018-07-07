import feastruct.post.results as results


class Node:
    """Class for a node to be used in finite element analyses.

    A node object is defined by its position in cartesian space and can store
    nodal displacement results from a static analysis, as well as eigenvector
    results from an eigenvalue analysis.

    :cvar int id: Unique node id
    :cvar coord: Cartesian coordinates of the node
    :vartype coord: list[float, float]
    :cvar dofs: Global degree of freedom numbers to be used in an analysis
    :vartype dofs: list[float, float, float]
    :cvar fixity: A list containing fixities for the nodal degrees of freedom
        for post-processing visualisation of nodal supports
    :vartype fixity: list[float, float, float]
    :cvar u: Nodal displacement results for static analyses
    :vartype u: :class:`feastruct.post.results.ResultList`
    :cvar buckling_v: Eigenvector results for buckling analyses
    :vartype buckling_v: :class:`feastruct.post.results.ResultList`
    :cvar frequency_v: Eigenvector results for frequency analyses
    :vartype frequency_v: :class:`feastruct.post.results.ResultList`
    """

    def __init__(self, id, coord):
        """Inits the Node class.

        :param int id: Unique node id
        :param coord: Cartesian coordinates of the node
        :type coord: list[float, float]
        """

        self.id = id
        self.coord = coord
        self.dofs = []
        self.fixity = [0, 0, 0]  # for post processing only
        self.u = results.ResultList()
        self.buckling_v = results.ResultList()
        self.frequency_v = results.ResultList()
        # TODO: check value types and int > 0

# TODO: document properties...

    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]

    @property
    def coords(self):
        return [self.x, self.y]

    def get_displacement(self, case_id):
        """Returns the displacement vector of the current node from the
        analysis defined by case_id.

        :param int case_id: Unique case id
        :return: Nodal displacement vector
        :rtype: :class:`numpy.ndarray`
        """

        return self.u.get_result(case_id).u

    def set_displacements(self, case_id, u):
        """Saves a displacement vector for the analysis defined by case_id.

        :param int case_id: Unique case id
        :param u: Nodal displacement vector
        :type u: :class:`numpy.ndarray`
        """

        self.u.set_result(results.Displacement(case_id, u))

    def get_buckling_results(self, case_id, buckling_mode):
        """Returns the eigenvalue (w) and nodal eigenvector (v) for the
        buckling analysis defined by case_id and the buckling mode defined by
        buckling_mode.

        :param int case_id: Unique case id
        :param int buckling_mode: Buckling mode number
        :return: (w, v)
        :rtype: tuple(float, :class:`numpy.ndarray`)
        """

        result = self.buckling_v.get_result(case_id, mode=buckling_mode)

        return (result.w, result.v)

    def set_buckling_results(self, case_id, buckling_mode, w, v):
        """Saves the buckling eigenvalue and nodal eigenvector for the buckling
        analysis defined by case_id and the buckling mode defined by
        buckling_mode.

        :param int case_id: Unique case id
        :param int buckling_mode: Buckling mode number
        :param float w: Buckling eigenvalue
        :param v: Buckling nodal eigenvector
        :type v: :class:`numpy.ndarray`
        """

        self.buckling_v.set_result(
            results.EigenResult(case_id, buckling_mode, w, v))

    def get_frequency_results(self, case_id, frequency_mode):
        """Returns the eigenvalue (w) and nodal eigenvector (v) for the
        frequency analysis defined by case_id and the frequency mode defined by
        frequency_mode.

        :param int case_id: Unique case id
        :param int frequency_mode: Frequency mode number
        :return: (w, v)
        :rtype: tuple(float, :class:`numpy.ndarray`)
        """

        result = self.frequency_v.get_result(case_id, mode=frequency_mode)

        return (result.w, result.v)

    def set_frequency_results(self, case_id, frequency_mode, w, v):
        """Saves the frequency eigenvalue and nodal eigenvector for the
        frequency analysis defined by case_id and the frequency mode defined by
        frequency_mode.

        :param int case_id: Unique case id
        :param int buckling_mode: Frequency mode number
        :param float w: Frequency eigenvalue
        :param v: Frequency nodal eigenvector
        :type v: :class:`numpy.ndarray`
        """

        self.frequency_v.set_result(
            results.EigenResult(case_id, frequency_mode, w, v))
