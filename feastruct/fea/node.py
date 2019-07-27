from operator import itemgetter


class Node:
    """Class for a node to be used in finite element analyses.

    A node object is defined by its position in 3D cartesian space and has six associated degrees
    of freedom.

    :cvar coords: Cartesian coordinates of the node
    :vartype coords: list[float, float, float]
    :cvar dofs: List of degrees of freedom for the node *(x, y, z, rx, ry, rz)*
    :vartype dofs: list[:class:`~feastruct.fea.node.DoF`]
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

    def get_dofs(self, node_dof_nums):
        """Returns the degree of freedom objects relating to the node_dof_nums list.

        :param node_dof_nums: Degree of freedom indices to return
        :type node_dof_nums: list[int]
        :returns: List containing degree of freedom objects
        :rtype: list[:class:`~feastruct.fea.node.DoF`]
        """

        # get all dofs corresponding with indices in 'node_dof_nums'
        dof_list = itemgetter(*node_dof_nums)(self.dofs)

        # ensure a list is returned
        return [dof_list] if type(dof_list) is not tuple else list(dof_list)


class DoF:
    """Class for a degree of freedom to be used in finite element analyses.

    A degree of freedom relates to a specific node and has a specific degree of freedom number used
    in the finite element analysis. Displacement results are stored within the degree of freedom.

    :cvar node: Parent node
    :vartype node: :class:`~feastruct.fea.node.Node`
    :cvar int node_dof_num: Node degree of freedom number
    :cvar int global_dof_num: Global degree of freedom number
    :cvar displacements: A list of displacement objects for the dof
    :vartype displacements: list[:class:`~feastruct.fea.node.Displacement`]
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

    def get_displacement(self, analysis_case):
        """Returns the displacement value relating to an :class:`~feastruct.fea.cases.AnalysisCase`.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :return: Displacement value
        :rtype: float

        :raises Exception: If a displacement corresponding to analysis_case cannot be found
        """

        for displacement in self.displacements:
            if analysis_case == displacement.analysis_case:
                return displacement.disp

        # if nothing is found
        str = 'Displacement corresponding to dof {0}'.format(self)
        str += ' could not be found for analysis case {0}'.format(analysis_case)
        raise Exception(str)


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
