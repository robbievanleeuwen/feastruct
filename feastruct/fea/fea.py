import numpy as np
from feastruct.fea.node import Node
import feastruct.fea.cases as Cases
from feastruct.post.results import ResultList, FrameForceVector
from feastruct.fea.exceptions import FEAInputError


class fea:
    """Parent class for a finite element analysis.

    Establishes a template for each different type of finite element analysis,
    e.g. frame, membrane etc. and provides a number of generic methods which
    are useful for all types of analyses.

    :cvar nodes: Nodes used in the finite element analysis.
    :vartype nodes: list[:class:`feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis.
    :vartype elements: list[:class:`feastruct.fea.fea.FiniteElement`]
    :cvar freedom_cases: Freedom cases used in the finite element analysis.
    :vartype freedom_cases: list[:class:`feastruct.fea.cases.FreedomCase`]
    :cvar load_cases: Load cases used in the finite element analysis.
    :vartype load_cases: list[:class:`feastruct.fea.cases.LoadCase`]
    :cvar analysis_cases: Analysis cases used in the finite element analysis.
    :vartype analysis_cases: list[:class:`feastruct.fea.cases.AnalysisCase`]
    :cvar bool non_linear: Boolean defining the type of analysis
    :cvar int dofs: Number of degrees of freedom per node for the current
        analysis type
    """

    def __init__(self, nodes, elements, freedom_cases, load_cases,
                 analysis_cases, non_linear, dofs):
        """Inits the fea class.

        :param nodes: List of nodes with which to initialise the class.
        :type nodes: list[:class:`feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class.
        :type elements: list[:class:`feastruct.fea.fea.FiniteElement`]
        :param freedom_cases: List of freedom cases with which to initialise
            the class.
        :type freedom_cases: list[:class:`feastruct.fea.cases.FreedomCase`]
        :param load_cases: List of load cases with which to initialise the
            class.
        :type load_cases: list[:class:`feastruct.fea.cases.LoadCase`]
        :param analysis_cases: List of analysis cases with which to initialise
            the class.
        :type analysis_cases: list[:class:`feastruct.fea.cases.AnalysisCase`]
        :param bool non_linear: Boolean defining the type of analysis
        :param int dofs: Number of degrees of freedom per node for the current
            analysis type
        """

        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes

        if elements is None:
            self.elements = []
        else:
            self.elements = elements

        if freedom_cases is None:
            self.freedom_cases = []
        else:
            self.freedom_cases = freedom_cases

        if load_cases is None:
            self.load_cases = []
        else:
            self.load_cases = load_cases

        if analysis_cases is None:
            self.analysis_cases = []
        else:
            self.analysis_cases = analysis_cases

        self.non_linear = non_linear
        self.dofs = dofs

    def add_node(self, id, coord):
        """Adds a node to the fea class.

        Adds a :class:`feastruct.fea.node.Node` object to the fea object.
        Checks for an exisiting id and ensures that a node is not already in
        the desired location.

        :param int id: Unique id identifying the node
        :param coord: Cartesian coordinates of the node
        :type coord: List[float, float]
        :raises FEAInputError: If the id already exists or if a node
            already exists at the desired location.

        The following example creates an analysis object and creates a node at
        x = 1 and y = 2::

            analysis = fea()
            analysis.add_node(id=1, coord=[1, 2])
        """

        # TODO: check that node id and location does not already exist
        self.nodes.append(Node(id, coord))

        # raise exception if duplicate added

    def add_element(self, element):
        """Adds an element to the fea class.

        Adds a :class:`feastruct.fea.fea.FiniteElement` object to the fea
        object. Checks for an exisiting id.

        :param element: Element to be added to the analysis object.
        :type element: :class:`feastruct.fea.fea.FiniteElement`
        :raises FEAInputError: If the specified element id already exists.
        """

        # TODO: check that element id does not already exist
        self.elements.append(element)

        # raise exception if duplicate added

    def add_freedom_case(self, id, items=None):
        """Adds a freedom case to the fea class.

        Adds a :class:`feastruct.fea.cases.FreedomCase` object to the fea
        object. Checks for an exisiting id.

        :param int id: Freedom case unique id
        :param items: List of items to initialise the freedom case with
        :type item: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
        :return: A freedom case object
        :rtype: :class:`feastruct.fea.cases.FreedomCase`
        :raises FEAInputError: If the freedom case id already exists.

        The following example creates an analysis object and creates a freedom
        case with and id of 1::

            analysis = fea()
            fc1 = analysis.add_freedom_case(id=1)
        """

        # TODO: check to see if id already exists and raise exception

        fc = Cases.FreedomCase(self, id, items)
        self.freedom_cases.append(fc)
        return fc

    def add_load_case(self, id, items=None):
        """Adds a load case to the fea class.

        Adds a :class:`feastruct.fea.cases.LoadCase` object to the fea object.
        Checks for an exisiting id.

        :param int id: Load case unique id
        :param items: List of items to initialise the load case with
        :type item: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
        :return: A load case object
        :rtype: :class:`feastruct.fea.cases.LoadCase`
        :raises FEAInputError: If the load case id already exists.

        The following example creates an analysis object and creates a load
        case with and id of 1::

            analysis = fea()
            lc1 = analysis.add_load_case(id=1)
        """

        # TODO: check to see if id already exists and raise exception

        lc = Cases.LoadCase(self, id, items)
        self.load_cases.append(lc)
        return lc

    def add_analysis_case(self, id, fc_id, lc_id):
        """Adds an analysis case to the fea class.

        Adds a :class:`feastruct.fea.cases.AnalysisCase` object to the fea
        object. Checks for an exisiting id.

        :param int id: Load case unique id
        :param int fc_id: Freedom case id to add to the analysis case
        :param int lc_id: Load case id to add to the analysis case
        :raises FEAInputError: If the load case id already exists.

        The following example creates an analysis object, a load and freedom
        case, and an analysis case using these newly created load and freedom
        cases::

            analysis = fea()
            fc1 = analysis.add_freedom_case(id=1)
            lc1 = analysis.add_load_case(id=1)
            analysis.add_analysis_case(id=1, fc_id=1, lc_id=1)
        """

        # TODO: check to see if id already exists and raise exception

        self.analysis_cases.append(Cases.AnalysisCase(self, id, fc_id, lc_id))

    def find_node(self, node_id):
        """Finds and returns the Node object given its node_id.

        :param int node_id: Unique node id
        :return: The Node object corresponding to the node_id
        :rtype: :class:`feastruct.fea.node.Node`
        :raises FEAInputError: If the Node corresponding to node_id cannot be
            found
        """

        # return the Node object whose id matches the given node_id
        for node in self.nodes:
            if node.id == node_id:
                return node
        else:
            raise FEAInputError("Cannot find node_id: {}".format(node_id))

    def get_node_lims(self):
        """Finds and returns the minimum and maximum x and y values within the
        current analysis.

        :return: (xmin, xmax, ymin, ymax)
        :rtype: tuple(float, float, float, float)
        """

        xmin = self.nodes[0].x
        xmax = self.nodes[0].x
        ymin = self.nodes[0].y
        ymax = self.nodes[0].y

        for node in self.nodes:
            xmin = min(xmin, node.x)
            xmax = max(xmax, node.x)
            ymin = min(ymin, node.y)
            ymax = max(ymax, node.y)

        return (xmin, xmax, ymin, ymax)

    def find_analysis_case(self, case_id):
        """Finds and returns the analysis case given its case_id

        :param int case_id: Unique case id
        :return: Analysis case corresponding to the case_id
        :rtype: :class:`feastruct.fea.cases.AnalysisCase`
        :raises FEAInputError: If the analysis case corresponding to case_id
            cannot be found
        """

        # return the AnalysisCase object whose id matches the given case_id
        for analysis_case in self.analysis_cases:
            if analysis_case.id == case_id:
                return analysis_case
        else:
            raise FEAInputError("Cannot find AnalysisCase id: {}".format(
                case_id))


class FiniteElement:
    """Parent class for a finite element.

    Establishes a template for each different type of finite element,
    e.g. frame element, membrane element, plate element etc. and provides a
    number of generic methods which are useful for all types of elements.

    :cvar analysis: Analysis object
    :vartype analysis: :class:`feastruct.fea.fea.fea`
    :cvar int id: Unique finite element id
    :cvar nodes: List of node objects that define the geometry of the finite
        element
    :vartype nodes: list[:class:`feastruct.fea.node.Node`]
    :cvar f_int: Internal force vector results for an arbitray number of
        analysis cases
    :vartype f_int: :class:`feastruct.post.results.ResultList`
    """

    def __init__(self, analysis, id, node_ids):
        """inits the FiniteElement class.

        :param analysis: Analysis object
        :type analysis: :class:`feastruct.fea.fea.fea`
        :param int id: Unique finite element id
        :param node_ids: A list of node ids defining the geometry of the
            element
        :type node_ids: list[int]
        """

        self.analysis = analysis
        self.id = id
        self.nodes = []
        self.f_int = ResultList()

        # add references to the node objects
        try:
            self.get_nodes(node_ids)
        except FEAInputError as error:
            print(error)

    def get_nodes(self, node_ids):
        """Finds and returns a list of node objects corresponding to nodes with
        ids in node_ids.

        :param node_ids: A list of unique node ids
        :type node_ids: list[int]
        :return: A list of node objects corresponding to the list of node_ids
        :rtype: list[:class:`feastruct.fea.node.Node`]
        :raises FEAInputError: If a node with node_id cannot be find in the
            analysis object
        """

        for node_id in node_ids:
            try:
                self.nodes.append(self.analysis.find_node(node_id))
            except FEAInputError as error:
                print("Error in FiniteElement id: {}. {}".format(
                    self.id, error))

    def get_node_coords(self):
        """Returns an array of the cartesian coordinates defining the geometry
        of the finite element.

        :return: An (n x 2) array of node coordinates, where n is the number of
            nodes for the given finite element
        :rtype: :class:`numpy.ndarray`
        """

        coord_list = []

        for node in self.nodes:
            coord_list.append(node.coords)

        return np.array(coord_list)

    def get_nodal_displacements(self, case_id):
        """Returns an array of the nodal displacements of each degree of
        freedom in the finite element for an analysis case defined by case_id.

        :param int case_id: Unique case id
        :return: An (n x ndof) array of degree of freedom displacements, where
            n is the number of nodes for the given finite element and ndof is
            the number of degrees of freedom per node for the current analysis
            type
        :rtype: :class:`numpy.ndarray`
        """

        disp_list = []  # allocate list of displacements

        for node in self.nodes:
            disp_list.append(node.get_displacement(case_id))

        return np.array(disp_list)

    def get_buckling_results(self, case_id, buckling_mode):
        """Returns the eigenvalue corresponding to the buckling analysis
        defined by the analysis case id and the buckling mode. Also returns an
        array of eigenvector values corresponding to each degree of freedom in
        the finite element for the given analysis case and buckling mode.

        :param int case_id: Unique case id
        :param int buckling_mode: Buckling mode number
        :return: (w, v) A tuple containing the eigenvalue 'w' and an (n x ndof)
            array of degree of freedom eigenvector values 'v', where n is the
            number of nodes for the given finite element and ndof is the number
            of degrees of freedom per node for the current analysis type
        :rtype: tuple(float, :class:`numpy.ndarray`)
        """

        v_list = []  # allocate list of eigenvectors

        for node in self.nodes:
            (w, v) = node.get_buckling_results(case_id, buckling_mode)
            v_list.append(v)

        return (w, np.array(v_list))

    def get_frequency_results(self, case_id, frequency_mode):
        """Returns the frequency corresponding to the natural frequency
        analysis defined by the analysis case id and the frequency mode. Also
        returns an array of eigenvector values corresponding to each degree of
        freedom in the finite element for the given analysis case and frequency
        mode.

        :param int case_id: Unique case id
        :param int frequency_mode: Frequency mode number
        :return: (w, v) A tuple containing the frequency 'w' and an (n x ndof)
            array of degree of freedom eigenvector values 'v', where n is the
            number of nodes for the given finite element and ndof is the number
            of degrees of freedom per node for the current analysis type
        :rtype: tuple(float, :class:`numpy.ndarray`)
        """

        v_list = []  # allocate list of eigenvectors

        for node in self.nodes:
            (w, v) = node.get_frequency_results(case_id, frequency_mode)
            v_list.append(v)

        return (w, np.array(v_list))

    def get_dofs(self):
        """Finds and returns an array of global degree of freedom numbers
        corresponding to the degrees of freedom of the finite element.

        :return: A integer array of degrees of freedom, with a length of
            n*ndof, where n is the number of nodes for the given finite element
            and ndof is the number of degrees of freedom per node for the
            current analysis type
        :rtype: :class:`numpy.ndarray`
        """

        dof_list = np.array([], dtype=np.int32)

        for node in self.nodes:
            dof_list = np.append(dof_list, node.dofs)

        return dof_list

    def get_fint(self, case_id):
        """Returns the internal force vector corresponding to the analysis case
        id.

        :param int case_id: Unique analysis case id
        :return: A result item containing the internal force vector
        :rtype: :class:`feastruct.post.results.ResultItem`
        """

        return self.f_int.get_result(case_id)

    def set_fint(self, case_id, f_int):
        """Saves the internal force vector result corresponding to the analysis
        case id.

        :param int case_id: Unique analysis case id
        :param f_int: Internal force vector array, with a length of n*ndof,
            where n is the number of nodes for the given finite element
            and ndof is the number of degrees of freedom per node for the
            current analysis type
        :type f_int: :class:`numpy.ndarray`
        """

        # save internal force vector in global coordinate system
        self.f_int.set_result(FrameForceVector(case_id, f_int))
