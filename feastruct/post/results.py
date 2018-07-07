from feastruct.fea.exceptions import FEAInputError


class ResultList:
    """Class for storing a list of result items.

    Stores a list of results in class variable self.results. Provides methods
    for getting and setting results for static and buckling/frequency analyses
    and overwrites results if the case id and/or mode match.

    :cvar results: List of result items
    :vartype results: list[:class:`feastruct.post.results.ResultItem`]
    """

    # TODO: check cvar results hyperlink in docs

    def __init__(self):
        """Inits the ResultList class.
        """

        self.results = []

    def get_result(self, case_id, mode=None):
        """Gets a result in the result list defined by the result case case_id
        and the mode for buckling/frequency analyses.

        :param int case_id: Unique case id
        :param int mode: Buckling/frequency mode number
        :return: Result defined by case_id (and) mode.
        :rtype: :class:`feastruct.post.results.ResultItem`
        :raises FEAInputError: If the case_id cannot be found in the list of
            results
        """

        # find case_id in results list
        for result in self.results:
            if result.case_id == case_id:
                if mode is None:
                    return result
                else:
                    if result.mode == mode:
                        return result

        # if case_id cannot be found
        raise FEAInputError("Cannot find a result for case_id: {}".format(
            case_id))

        # TODO: raise error if mode number cannot be found

    def set_result(self, new_result):
        """Sets a result in the result list. Replaces the result if the case id
        already exists and only replaces results for buckling/frequency
        analyses if both the case id and mode already exist. Appends the result
        if the case id does not yet exist.

        :param new_result: Result to be appended to the result list
        :type new_result: :class:`feastruct.post.results.ResultItem`
        """

        # check to see if case_id already exists
        for result in self.results:
            if result.case_id == new_result.case_id:
                # if the ResultItem has a mode (i.e. is an EigenResult)
                if type(result) is EigenResult:
                    # replace result if the mode matches
                    if result.mode == new_result.mode:
                        result = new_result
                        return
                    else:
                        # append to result list if this is a new mode
                        self.results.append(new_result)

                else:
                    # replace result if case_id already exists
                    result = new_result
                    return

        # append to result list if this is a new case_id
        self.results.append(new_result)


class ResultItem:
    """Parent class for a finite element result.

    :cvar int case_id: Unique case id
    """

    def __init__(self, case_id):
        """Inits the ResultItem class.

        :param int case_id: Unique case id
        """

        self.case_id = case_id


class Displacement(ResultItem):
    """Class for a displacement result.

    :cvar int case_id: Unique case id
    :cvar u: Nodal displacement vector with length equal to the number of
        degrees of freedom
    :vartype u: list[float]
    """

    def __init__(self, case_id, u):
        """Inits the Displacement class.

        :param int case_id: Unique case id
        :param u: Nodal displacement vector with length equal to the number of
            degrees of freedom
        :vartype u: list[float]
        """

        super().__init__(case_id)
        self.u = u


class EigenResult(ResultItem):
    """Class for an eigenvalue and eigenvector result.

    :cvar int case_id: Unique case id
    :cvar int mode: Buckling/frequency mode number
    :cvar float w: Buckling load factor or natural frequency
    :cvar v: Nodal eigenvector with length equal to the number of degrees of
        freedom
    :vartype v: list[float]
    """

    def __init__(self, case_id, mode, w, v):
        """Inits the EigenResult class.

        :param int case_id: Unique case id
        :param int mode: Buckling/frequency mode number
        :param float w: Buckling load factor or natural frequency
        :param v: Nodal eigenvector with length equal to the number of degrees
            of freedom
        :type v: list[float]
        """

        super().__init__(case_id)
        self.mode = mode
        self.w = w
        self.v = v


class Force(ResultItem):
    """Class for a scalar force result e.g. reaction force.

    :cvar int case_id: Unique case id
    :cvar float f: Scalar force result
    """

    def __init__(self, case_id, f):
        """Inits the Force class.

        :param int case_id: Unique case id
        :param float f: Scalar force result
        """

        super().__init__(case_id)
        self.f = f


class FrameForceVector(ResultItem):
    """Class for a frame element internal force vector. The sign conventions
    are as follows:

    * Axial Force: Tension positive
    * Shear Force: Clockwise shear deformation positive
    * Bending Moment: Clockwise bending moment positive

    :cvar int case_id: Unique case id
    :cvar float N1: Internal axial force at node 1
    :cvar float V1: Internal shear force at node 1
    :cvar float M1: Internal bending moment at node 1
    :cvar float N2: Internal axial force at node 2
    :cvar float V2: Internal shear force at node 2
    :cvar float M2: Internal bending moment at node 2
    """

    def __init__(self, case_id, f):
        """Inits the FrameForceVector class.

        :param int case_id: Unique case id
        :param f: Frame internal force vector with length (n x ndof) where n is
            the number of nodes for the given finite element and ndof is the
            number of degrees of freedom per node for the current analysis type
        :type f: list[float]
        """

        super().__init__(case_id)
        self.N1 = f[0]
        self.V1 = f[1]
        self.M1 = f[2]
        self.N2 = f[3]
        self.V2 = f[4]
        self.M2 = f[5]
