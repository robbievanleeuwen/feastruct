import numpy as np
from feastruct.fea.node import Node


class FiniteElementAnalysis:
    """Parent class for a finite element analysis.

    Establishes a template for each different type of finite element analysis, e.g. frame,
    membrane, plate etc. and provides a number of generic methods which are useful for all types of
    analyses.

    :cvar nodes: Nodes used in the finite element analysis
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis
    :vartype elements: list[:class:`~feastruct.fea.fea.FiniteElement`]
    :cvar nfa: Node freedom arrangement
    :vartype nfa: list[bool]
    """

    def __init__(self, nodes, elements, nfa):
        """Inits the FiniteElementAnalysis class.

        :param nodes: List of nodes with which to initialise the class
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class
        :type elements: list[:class:`~feastruct.fea.fea.FiniteElement`]
        :param nfa: Node freedom arrangement
        :type nfa: list[bool]
        """

        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes

        if elements is None:
            self.elements = []
        else:
            self.elements = elements

        self.nfa = nfa

    def create_node(self, coords):
        """Creates a node and adds it to the Fea object.

        Creates and returns a :class:`~feastruct.fea.node.Node` object and adds it to the
        :class:`~feastruct.fea.fea.Fea` object.

        :param coords: Cartesian coordinates of the node *([x], [x, y] or [x, y, z])*
        :type coords: list[float]
        :returns: Node object
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
        :returns: Element object
        :rtype: :class:`~feastruct.fea.fea.FiniteElement`
        """

        # TODO:catch any errors thrown by element creation

        self.elements.append(element)

        return element

    def get_node_lims(self):
        """Finds and returns the minimum and maximum x, y and z values within the current analysis.

        :returns: (xmin, xmax, ymin, ymax, zmin, zmax)
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
    :cvar efs: Element freedom signature
    :vartype efs: list[bool]
    :cvar f_int: List of internal force vector results stored for each analysis case
    :vartype f_int: list[:class:`~feastruct.fea.fea.ForceVector`]
    """

    def __init__(self, nodes, material, efs):
        """Inits the FiniteElement class.

        :param nodes: List of node objects that define the geometry of the finite element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param efs: Element freedom signature
        :type efs: list[bool]
        """

        self.nodes = nodes
        self.material = material
        self.efs = efs
        self.f_int = []

    def get_node_coords(self):
        """Returns a NumPy array of the cartesian coordinates defining the geometry of the finite
        element.

        :returns: An *(n x 3)* array of node coordinates, where *n* is the number of nodes for the
            given finite element
        :rtype: :class:`numpy.ndarray`
        """

        coord_list = []

        # loop through all nodes of the element
        for node in self.nodes:
            # append the node coordinates to the list
            coord_list.append(node.coords)

        return np.array(coord_list)

    def get_ndof(self):
        """Returns the number of active degrees of freedom for the element.

        :returns: Number of active degrees of freedom for the element
        :rtype: int
        """

        return len(self.nodes) * sum(self.efs)

    def get_dofs(self):
        """Finds and returns a list of DoF objects corresponding to the degrees of freedom of the
        finite element.

        :returns: A list of DoF objects, with a length of *(n_nodes x n_dof)*
        :rtype: list[list[:class:`~feastruct.fea.node.DoF`]]
        """

        dof_list = []

        # loop through all nodes of the element
        for node in self.nodes:
            # append the node dofs to the list
            dof_list.append(node.get_dofs(freedom_signature=self.efs))

        return dof_list

    def get_gdof_nums(self):
        """Returns an array of global degree of freedom numbers corresponding to the degrees of
        freedom of the finite element.

        :returns: A integer array of global degrees of freedom, with a length of *(1 x n)*, where n
            is *n_nodes x n_dof*
        :rtype: :class:`numpy.ndarray`
        """

        # get a list of dof objects for all nodes
        dof_list = self.get_dofs()

        # allocate a list of ints for the global
        gdof_num_list = []

        # loop through nodes in the dof list
        for node in dof_list:
            # loop through dofs in the current node
            for dof in node:
                gdof_num_list.append(dof.global_dof_num)

        return np.array(gdof_num_list, dtype=np.int32)

    def apply_nfa(self):
        """Applies the element freedom signature to all nodes in the element to generate a node
        freedom allocation.
        """

        # loop through all the nodes in the element
        for node in self.nodes:
            # loop through all the freedoms in the current node signature
            for (i, nf) in enumerate(node.nfs):
                # apply the element freedom signature
                node.nfs[i] = nf or self.efs[i]

    def get_nodal_displacements(self, analysis_case):
        """Returns an array of the nodal displacements for each degree of freedom in the finite
        element for the analysis_case.

        :param analysis_case: Analysis case relating to the displacement
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :returns: An *(n_nodes x n_dof)* array of degree of freedom displacements
        :rtype: :class:`numpy.ndarray`
        """

        disp_list = []  # allocate list of displacements

        # get all the dof objects for the element
        dof_list = self.get_dofs()

        # loop through each node's dofs
        for node_dofs in dof_list:
            # list of displacements for the current node
            node_list = []

            # loop through each dof in the node
            for dof in node_dofs:
                # add the dof displacement to the node list
                node_list.append(dof.get_displacement(analysis_case))

            disp_list.append(node_list)

        return np.array(disp_list)

    def get_buckling_results(self, analysis_case, buckling_mode=0):
        """Returns the eigenvalue corresponding to the buckling analysis defined by the
        analysis_case and the buckling_mode. Also returns an array of eigenvector values
        corresponding to each degree of freedom in the finite element.

        :param analysis_case: Analysis case relating to the buckling analysis
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int buckling_mode: Buckling mode number

        :returns: Eigenvalue and eigenvectors *(w, v)*. A tuple containing the eigenvalue *w* and
            an *(n_nodes x n_dof)* array of eigenvector values *v*.
        :rtype: tuple(float, :class:`numpy.ndarray`)
        """

        v_list = []  # allocate list of eigenvectors

        # get all the dof objects for the element
        dof_list = self.get_dofs()

        # loop through each node's dofs
        for node_dofs in dof_list:
            # list of eigenvectors for the current node
            node_list = []

            # loop through each dof in the node
            for dof in node_dofs:
                (w, v) = dof.get_buckling_mode(
                    analysis_case=analysis_case, buckling_mode=buckling_mode)

                # add the dof eigenvector to the node list
                node_list.append(v)

            v_list.append(node_list)

        return (w, np.array(v_list))

    def get_frequency_results(self, analysis_case, frequency_mode=0):
        """Returns the eigenvalue corresponding to the frequency analysis defined by the
        analysis_case and the frequency_mode. Also returns an array of eigenvector values
        corresponding to each degree of freedom in the finite element.

        :param analysis_case: Analysis case relating to the frequency analysis
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int frequency_mode: Frequency mode number

        :returns: Eigenvalue and eigenvectors *(w, v)*. A tuple containing the eigenvalue *w* and
            an *(n_nodes x n_dof)* array of eigenvector values *v*.
        :rtype: tuple(float, :class:`numpy.ndarray`)
        """

        v_list = []  # allocate list of eigenvectors

        # get all the dof objects for the element
        dof_list = self.get_dofs()

        # loop through each node's dofs
        for node_dofs in dof_list:
            # list of eigenvectors for the current node
            node_list = []

            # loop through each dof in the node
            for dof in node_dofs:
                (w, v) = dof.get_frequency_mode(
                    analysis_case=analysis_case, frequency_mode=frequency_mode)

                # add the dof eigenvector to the node list
                node_list.append(v)

            v_list.append(node_list)

        return (w, np.array(v_list))

    def get_fint(self, analysis_case):
        """Returns the internal force vector relating to an
        :class:`~feastruct.fea.cases.AnalysisCase`.

        :param analysis_case: Analysis case relating to the internal force vector
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :returns: Internal force vector
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

    def get_shape_function(self, xi):
        """Placeholder for the get_shape_function method.

        Returns the value of the shape functions at *xi*.

        :param float xi: Position along the element

        :returns: Value of the shape functions at *xi*
        :rtype: :class:`numpy.ndarray`
        """

        pass

    def get_stiffness_matrix(self):
        """Placeholder for the get_stiffness_matrix method.

        Gets the stiffness matrix for a FiniteElement.

        :returns: 6 x 6 element stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        pass

    def get_geometric_stiff_matrix(self, analysis_case):
        """Placeholder for the get_geometric_stiff_matrix method.

        Gets the geometric stiffness matrix for a FiniteElement.

        :param analysis_case: Analysis case from which to extract the axial force
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: Element geometric stiffness matrix
        :rtype: :class:`numpy.ndarray`
        """

        pass

    def get_mass_matrix(self):
        """Placeholder for the get_mass_matrix method.

        Gets the mass matrix for a for a FiniteElement.

        :returns: Element mass matrix
        :rtype: :class:`numpy.ndarray`
        """

        pass

    def get_internal_actions(self, analysis_case):
        """Returns the internal actions for a FiniteElement.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`

        :returns: An array containing the internal actions for the element
        :rtype: :class:`numpy.ndarray`
        """

        pass


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
