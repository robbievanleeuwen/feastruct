class Node:
    """Class for a node to be used in finite element analyses.

    A node object is defined by its position in 3D cartesian space and has six associated degrees
    of freedom.

    :cvar coords: Cartesian coordinates of the node
    :vartype coords: list[float, float, float]
    :cvar dofs: List of degrees of freedom for the node *(x, y, z, rx, ry, rz)*
    :vartype dofs: list[:class:`~feastruct.fea.node.DoF`]
    :cvar nfs: Node freedom signature
    :vartype nfs: list[bool]
    """

    def __init__(self, coords):
        """Inits the Node class.

        :param coords: Cartesian coordinates of the node *([x], [x, y] or [x, y, z])*
        :type coords: list[float]
        """

        if len(coords) == 1:
            self.coords = [coords[0], 0, 0]
        elif len(coords) == 2:
            self.coords = [coords[0], coords[1], 0]
        elif len(coords) == 3:
            self.coords = coords
        else:
            # TODO: throw an error
            return

        # create six dofs
        self.dofs = []

        for i in range(6):
            self.dofs.append(DoF(node=self, node_dof_num=i))

        # set the node freedom signature
        self.nfs = [False, False, False, False, False, False]

    # get x-coordinate
    @property
    def x(self):
        return self.coords[0]

    # get y-coordinate
    @property
    def y(self):
        return self.coords[1]

    # get z-coordinate
    @property
    def z(self):
        return self.coords[2]

    def move_node(self, vec):
        """Moves the current node by the vector vec *([dx], [dx, dy] or [dx, dy, dz])*.

        :param vec: Vector to move the current node by
        :type vec: list[float]
        """

        if len(vec) == 1:
            self.coords[0] += vec[0]
        elif len(vec) == 2:
            self.coords[0] += vec[0]
            self.coords[1] += vec[1]
        elif len(vec) == 3:
            self.coords += vec
        else:
            # TODO: throw an error
            return

    def copy_node(self, vec):
        """Copies and returns a new node shifted by the vector vec
        *([dx], [dx, dy] or [dx, dy, dz])*.

        :param vec: Vector to move the current node by
        :type vec: list[float]

        :returns: New node object
        :rtype: :class:`~feastruct.fea.node.Node`
        """

        new_coords = self.coords

        if len(vec) == 1:
            new_coords[0] += vec[0]
        elif len(vec) == 2:
            new_coords[0] += vec[0]
            new_coords[1] += vec[1]
        elif len(vec) == 3:
            new_coords += vec
        else:
            # TODO: throw an error
            return

        return Node(coords=new_coords)

    def get_dofs(self, freedom_signature):
        """Returns the degree of freedom objects relating to the freedom_signature list.

        :param freedom_signature: Degree of freedom indices to return
        :type freedom_signature: list[bool]
        :returns: List containing degree of freedom objects
        :rtype: list[:class:`~feastruct.fea.node.DoF`]
        """

        # initialise dof_list
        dof_list = []

        # loop through freedom_signature list
        for (i, freedom) in enumerate(freedom_signature):
            # if the freedom is True
            if freedom:
                dof_list.append(self.dofs[i])

        return dof_list

    def get_displacements(self, analysis_case):
        """Returns the displacement vector for the current node and analysis case. If the degree
        of freedom is not used in the current analysis case, a *None* value is assigned.

        Degree of freedoms are ordered as follows: *(x, y, z, rx, ry, rz)*

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Displacement vector of length *(1 x 6)*
        :rtype: float
        """

        # assign displacement vector
        disp_vector = []

        # loop throught the degrees of freedom
        for dof in self.dofs:
            try:
                disp = dof.get_displacement(analysis_case=analysis_case)
            except Exception:
                disp = None

            disp_vector.append(disp)

        return disp_vector


class DoF:
    """Class for a degree of freedom to be used in finite element analyses.

    A degree of freedom relates to a specific node and has a specific degree of freedom number used
    in the finite element analysis. Displacement results are stored within the degree of freedom.

    :cvar node: Parent node
    :vartype node: :class:`~feastruct.fea.node.Node`
    :cvar int node_dof_num: Node degree of freedom number
    :cvar int global_dof_num: Global degree of freedom number
    :cvar displacements: A list of Displacement objects for the dof
    :vartype displacements: list[:class:`~feastruct.fea.node.Displacement`]
    :cvar buckling_results: A list of BucklingModes objects for the dof
    :vartype buckling_results: list[:class:`~feastruct.fea.node.BucklingMode`]
    :cvar frequency_results: A list of FrequencyModes objects for the dof
    :vartype frequency_results: list[:class:`~feastruct.fea.node.FrequencyMode`]
    """

    def __init__(self, node, node_dof_num):
        """Inits the DoF class.

        :param node: Parent node
        :type node: :class:`~feastruct.fea.node.Node`
        """

        self.node = node
        self.node_dof_num = node_dof_num
        self.global_dof_num = None
        self.displacements = []
        self.buckling_results = []
        self.frequency_results = []

    def get_displacement(self, analysis_case):
        """Returns the displacement value relating to an :class:`~feastruct.fea.cases.AnalysisCase`.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Displacement value
        :rtype: float

        :raises Exception: If a displacement corresponding to analysis_case cannot be found
        """

        # loop through all the displacements
        for displacement in self.displacements:
            # if the analysis_case is found
            if analysis_case == displacement.analysis_case:
                # return the displacement
                return displacement.disp

        # if nothing is found raise an exception
        str = 'Displacement corresponding to dof {0}'.format(self)
        str += ' could not be found for analysis case {0}'.format(analysis_case)
        raise Exception(str)

    def get_buckling_mode(self, analysis_case, buckling_mode):
        """Returns the value of the eigenvalue and eigenvector for a given buckling_mode related to
        an :class:`~feastruct.fea.cases.AnalysisCase`.

        :param analysis_case: Analysis case relating to the buckling results
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int buckling_mode: Buckling mode number

        :returns: Eigenvalue and eigenvector corresponding to an analysis_case and buckling_mode
        :rtype: tuple(float, float)

        :raises Exception: If buckling results corresponding to an analysis_case cannot be found,
            or if the buckling mode does not exist
        """

        # loop through all the buckling results
        for buckling_result in self.buckling_results:
            # if the analysis case is found
            if analysis_case == buckling_result.analysis_case:
                # loop through all the modes
                for (i, mode) in enumerate(buckling_result.buckling_modes):
                    # if the mode is found
                    if mode == buckling_mode:
                        # return the eigenvalue and eigevector
                        return (buckling_result.w[i], buckling_result.v[i])
                # if the mode is not found
                else:
                    str = 'Buckling result for dof {0} corresponding to buckling mode {1}'.format(
                        self, buckling_mode)
                    str += ' cannot be found for analysis case {0}.'.format(analysis_case)
                    str += ' Analysis case located, but buckling mode not found.'
                    raise Exception(str)

        # if the analysis case is not found
        str = 'Buckling result for dof {0} could not be found for analysis case {1}'.format(
            self, analysis_case)
        raise Exception(str)

    def get_frequency_mode(self, analysis_case, frequency_mode):
        """Returns the value of the eigenvalue and eigenvector for a given frequency_mode related
        to an :class:`~feastruct.fea.cases.AnalysisCase`.

        :param analysis_case: Analysis case relating to the frequency results
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int frequency_mode: Frequency mode number

        :returns: Eigenvalue and eigenvector corresponding to an analysis_case and frequency_mode
        :rtype: tuple(float, float)

        :raises Exception: If frequency results corresponding to an analysis_case cannot be found,
            or if the frequency mode does not exist
        """

        # loop through all the frequency results
        for frequency_result in self.frequency_results:
            # if the analysis case is found
            if analysis_case == frequency_result.analysis_case:
                # loop through all the modes
                for (i, mode) in enumerate(frequency_result.frequency_modes):
                    # if the mode is found
                    if mode == frequency_mode:
                        # return the eigenvalue and eigevector
                        return (frequency_result.w[i], frequency_result.v[i])
                # if the mode is not found
                else:
                    str = 'Frequency result for dof {0} corresponding to frequency mode'.format(
                        self)
                    str += ' {1} cannot be found for analysis case {0}.'.format(
                        analysis_case, frequency_mode)
                    str += ' Analysis case located, but frequency mode not found.'
                    raise Exception(str)

        # if the analysis case is not found
        str = 'Frequency result for dof {0} could not be found for analysis case {1}'.format(
            self, analysis_case)
        raise Exception(str)

    def save_displacement(self, disp, analysis_case):
        """Saves a displacement result to the DoF.

        :param float disp: Value of the displacement
        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # check to see if there is already a displacement for the current analysis_case
        for displacement in self.displacements:
            if displacement.analysis_case == analysis_case:
                displacement.disp = disp
                return

        self.displacements.append(Displacement(disp=disp, analysis_case=analysis_case))

    def save_buckling_modes(self, buckling_modes, w, v, analysis_case):
        """Saves buckling results to the DoF.

        :param buckling_modes: Buckling mode numbers
        :type buckling_modes: list[int]
        :param w: Eigenvalues corresponding to the given modes
        :type w: list[float]
        :param v: Value of the eigevectors at the DoF corresponding to the given modes
        :type v: list[float]
        :param analysis_case: Analysis case relating to the buckling analysis
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # create buckling mode result
        new_result = BucklingModes(
            buckling_modes=buckling_modes, w=w, v=v, analysis_case=analysis_case
        )

        # loop through all the buckling results
        for buckling_result in self.buckling_results:
            # if the analysis case already exists, overwrite the results
            if buckling_result.analysis_case == analysis_case:
                buckling_result = new_result
                return

        # otherwise add a new result
        self.buckling_results.append(new_result)

    def save_frequency_modes(self, frequency_modes, w, v, analysis_case):
        """Saves frequency results to the DoF.

        :param frequency_modes: Frequency mode numbers
        :type frequency_modes: list[int]
        :param w: Eigenvalues corresponding to the given modes
        :type w: list[float]
        :param v: Value of the eigevectors at the DoF corresponding to the given modes
        :type v: list[float]
        :param analysis_case: Analysis case relating to the frequency analysis
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # create frequency mode result
        new_result = FrequencyModes(
            frequency_modes=frequency_modes, w=w, v=v, analysis_case=analysis_case
        )

        # loop through all the frequency results
        for frequency_result in self.frequency_results:
            # if the analysis case already exists, overwrite the results
            if frequency_result.analysis_case == analysis_case:
                frequency_result = new_result
                return

        # otherwise add a new result
        self.frequency_results.append(new_result)


class Displacement:
    """Class for storing a displacement at a degree of freedom for a specific analysis case.

    :cvar float disp: Displacement at a degree of freedom
    :cvar analysis_case: Analysis case relating to the displacement
    :vartype analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    """

    def __init__(self, disp, analysis_case):
        """Inits the Displacement class.

        :param float disp: Displacement at a degree of freedom
        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        self.disp = disp
        self.analysis_case = analysis_case


class BucklingModes:
    """Class for storing buckling modes at a degree of freedom for a specific analysis case.

    :cvar buckling_modes: Buckling mode numbers
    :vartype buckling_modes: list[int]
    :cvar w: Eigenvalues corresponding to the given modes
    :vartype w: list[float]
    :cvar float v: Value of the eigevectors at the DoF corresponding to the given modes
    :vartype v: list[float]
    :cvar analysis_case: Analysis case relating to the buckling analysis
    :vartype analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    """

    def __init__(self, buckling_modes, w, v, analysis_case):
        """Inits the Displacement class.

        :param buckling_modes: Buckling mode numbers
        :type buckling_modes: list[int]
        :param w: Eigenvalues corresponding to the given modes
        :type w: list[float]
        :param v: Value of the eigevectors at the DoF corresponding to the given modes
        :type v: list[float]
        :param analysis_case: Analysis case relating to the buckling analysis
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        self.buckling_modes = buckling_modes
        self.w = w
        self.v = v
        self.analysis_case = analysis_case


class FrequencyModes:
    """Class for storing frequency modes at a degree of freedom for a specific analysis case.

    :cvar frequency_modes: Frequency mode numbers
    :vartype frequency_modes: list[int]
    :cvar w: Eigenvalues corresponding to the given modes
    :vartype w: list[float]
    :cvar float v: Value of the eigevectors at the DoF corresponding to the given modes
    :vartype v: list[float]
    :cvar analysis_case: Analysis case relating to the frequency analysis
    :vartype analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    """

    def __init__(self, frequency_modes, w, v, analysis_case):
        """Inits the Displacement class.

        :param frequency_modes: Frequency mode numbers
        :type frequency_modes: list[int]
        :param w: Eigenvalues corresponding to the given modes
        :type w: list[float]
        :param v: Value of the eigevectors at the DoF corresponding to the given modes
        :type v: list[float]
        :param analysis_case: Analysis case relating to the frequency analysis
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        self.frequency_modes = frequency_modes
        self.w = w
        self.v = v
        self.analysis_case = analysis_case
