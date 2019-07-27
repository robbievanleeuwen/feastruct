import numpy as np
from feastruct.fea.node import Node
from feastruct.post.post import PostProcessor


class FiniteElementAnalysis:
    """Parent class for a finite element analysis.

    Establishes a template for each different type of finite element analysis, e.g. frame,
    membrane, plate etc. and provides a number of generic methods which are useful for all types of
    analyses.

    :cvar nodes: Nodes used in the finite element analysis
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis
    :vartype elements: list[:class:`~feastruct.fea.fea.FiniteElement`]
    :cvar int dims: Number of dimensions used for the current analysis type
    :cvar dofs: List of the degrees of freedom used in the current analysis type
    :vartype dofs: list[int]
    :cvar post: Post-processor object
    :vartype post: :class:`feastruct.post.post.PostProcessor`
    """

    def __init__(self, nodes, elements, dims, dofs):
        """Inits the FiniteElementAnalysis class.

        :param nodes: List of nodes with which to initialise the class
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class
        :type elements: list[:class:`~feastruct.fea.fea.FiniteElement`]
        :param int dims: Number of dimensions used for the current analysis type
        :param dofs: List of the degrees of freedom used in the current analysis type
        :type dofs: list[int]
        """

        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes

        if elements is None:
            self.elements = []
        else:
            self.elements = elements

        self.dims = dims
        self.dofs = dofs
        self.post = PostProcessor(self)

    def create_node(self, coords):
        """Creates a node and adds it to the Fea object.

        Creates and returns a :class:`~feastruct.fea.node.Node` object and adds it to the
        :class:`~feastruct.fea.fea.Fea` object.

        :param coords: Cartesian coordinates of the node *([x], [x, y] or [x, y, z])*
        :type coords: list[float]
        :return: Node object
        :rtype: :class:`~feastruct.fea.node.Node`
        """

        # TODO:catch any errors thrown by node creation

        new_node = Node(coords)
        self.nodes.append(new_node)

        return new_node

    def create_element(self, element):
        """Creates a finite element and adds it to the Fea object.

        Creates and returns a :class:`~feastruct.fea.fea.FiniteElement` object ands adds it to the
        :class:`~feastruct.fea.fea.Fea` object.

        :param element: Element to be added to the analysis object
        :type element: :class:`~feastruct.fea.fea.FiniteElement`
        :return: Element object
        :rtype: :class:`~feastruct.fea.fea.FiniteElement`
        """

        # TODO:catch any errors thrown by element creation

        self.elements.append(element)

        return element

    def get_node_lims(self):
        """Finds and returns the minimum and maximum x, y and z values within the current analysis.

        :return: (xmin, xmax, ymin, ymax, zmin, zmax)
        :rtype: tuple(float)
        """

        for (i, node) in enumerate(self.nodes):
            if i == 0:
                xmin = node.x
                xmax = node.x
                ymin = node.y
                ymax = node.y
                zmin = node.z
                zmax = node.z

            xmin = min(xmin, node.x)
            xmax = max(xmax, node.x)
            ymin = min(ymin, node.y)
            ymax = max(ymax, node.y)
            zmin = min(zmin, node.z)
            zmax = max(zmax, node.z)

        return (xmin, xmax, ymin, ymax, zmin, zmax)


class FiniteElement:
    """Parent class for a finite element.

    Establishes a template for each different type of finite element, e.g. frame element, membrane
    element, plate element etc. and provides a number of generic methods which are useful for all
    types of elements.

    :cvar nodes: List of node objects that define the geometry of the finite element
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar material: Material object for the element
    :vartype material: :class:`~feastruct.pre.material.Material`
    :cvar f_int: List of internal force vector results stored for each analysis case
    :vartype f_int: list[:class:`~feastruct.fea.fea.ForceVector`]
    """

    def __init__(self, nodes, material):
        """Inits the FiniteElement class.

        :param nodes: List of node objects that define the geometry of the finite element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        """

        self.nodes = nodes
        self.material = material
        self.f_int = []

    def get_node_coords(self):
        """Returns a NumPy array of the cartesian coordinates defining the geometry of the finite
        element.

        :return: An *(n x 3)* array of node coordinates, where *n* is the number of nodes for the
            given finite element
        :rtype: :class:`numpy.ndarray`
        """

        coord_list = []

        # loop through all nodes of the element
        for node in self.nodes:
            # append the node coordinates to the list
            coord_list.append(node.coords)

        return np.array(coord_list)

    def get_dofs(self, dof_nums):
        """Finds and returns a list of DoF objects corresponding to the degrees of freedom of the
        finite element.

        :param dof_nums: List of degrees of freedom used by the current finite element type e.g.
            [0, 1] for *x* and *y* translation
        :return: A list of DoF objects, with a length of *(n_nodes x n_dof)*
        :rtype: list[list[:class:`~feastruct.fea.node.DoF`]]
        """

        dof_list = []

        # loop through all nodes of the element
        for node in self.nodes:
            # append the node dofs to the list
            dof_list.append(node.get_dofs(dof_nums))

        return dof_list

    def get_gdof_nums(self, dof_nums):
        """Returns an array of global degree of freedom numbers corresponding to the degrees of
        freedom of the finite element.

        :param dof_nums: List of degrees of freedom used by the current finite element type e.g.
            [0, 1] for *x* and *y* translation
        :return: A integer array of global degrees of freedom, with a length of *(1 x n)*, where n
            is *n_nodes x n_dof*
        :rtype: :class:`numpy.ndarray`
        """

        # get a list of dof objects for all nodes
        dof_list = self.get_dofs(dof_nums)

        # allocate a list of ints for the global
        gdof_num_list = []

        # loop through nodes in the dof list
        for node in dof_list:
            # loop through dofs in the current node
            for dof in node:
                gdof_num_list.append(dof.global_dof_num)

        return np.array(gdof_num_list, dtype=np.int32)

    def get_nodal_displacements(self, dof_nums, analysis_case):
        """Returns an array of the nodal displacements for each degree of freedom in the finite
        element for the analysis_case.

        :param dof_nums: List of degrees of freedom used by the current finite element type e.g.
            [0, 1] for *x* and *y* translation
        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :return: An *(n_nodes x n_dof)* array of degree of freedom displacements
        :rtype: :class:`numpy.ndarray`
        """

        disp_list = []  # allocate list of displacements
        dof_list = self.get_dofs(dof_nums)

        for node_dofs in dof_list:
            node_list = []

            for dof in node_dofs:
                node_list.append(dof.get_displacement(analysis_case))

            disp_list.append(node_list)

        return np.array(disp_list)

    def get_fint(self, analysis_case):
        """Returns the internal force vector relating to an
        :class:`~feastruct.fea.cases.AnalysisCase`.

        :param analysis_case: Analysis case relating to the internal force vector
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :return: Internal force vector
        :rtype: :class:`numpy.ndarray`

        :raises Exception: If the force vector cannot be found for the analysis_case
        """

        for f_int in self.f_int:
            if analysis_case == f_int.analysis_case:
                return f_int.f

        # if nothing is found
        str = 'Force vector corresponding to element {0}'.format(self)
        str += ' could not be found for analysis case {0}'.format(analysis_case)
        raise Exception(str)

    def save_fint(self, f, analysis_case):
        """Adds an internal force vector result to the :class:`~feastruct.fea.fea.FiniteElement`.

        :param f: Internal force vector
        :type f: :class:`numpy.ndarray`
        :param analysis_case: Analysis case relating to the internal force vector
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        self.f_int.append(ForceVector(f, analysis_case))


class ForceVector:
    """Class for storing a force vector for a finite element.

    :cvar f: Force vector *(n x 1)*, where *n* is the number of dofs within the element
    :vartype: :class:`numpy.ndarray`
    :cvar analysis_case: Analysis case relating to the displacement
    :vartype analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    """

    def __init__(self, f, analysis_case):
        """Inits the ForceVector class.

        :cvar f: Force vector *(n x 1)*, where *n* is the number of dofs within the element
        :vartype: :class:`numpy.ndarray`
        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        self.f = f
        self.analysis_case = analysis_case
