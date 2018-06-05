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

    Attributes:
        nodes:          list of 'Node' items with which to initialise the
                        class.
        elements:       list of 'FiniteElement' items with which to initialise
                        the class.
        freedom_cases:  list of 'FreedomCase' itmes with which to initialise
                        the class.
        load_cases:     list of 'LoadCase' itmes with which to initialise
                        the class.
        analysis_cases: list of 'AnalysisCase' items with which to initialise
                        the class.
        non_linear:     boolean defining the type of analysis
                        (geometric & material).
        dofs:           structural degrees of freedom per node for current
                        analysis type.
    """

    def __init__(self, nodes, elements, freedom_cases, load_cases,
                 analysis_cases, non_linear, dofs):
        """Inits the fea class with..."""

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

        Adds a Node object to the fea object. Checks for an exisiting id and
        ensures that a node is not already in the desired location.

        Args:
            id:     unique integer identifying the node
            coord:  cartesian coordinates of the node i.e. [x, y]

        Returns:
            void

        Raises:
            FEAInputError:  Raised if the id already exists or if a node
                            already exists at the desired location.
        """

        # TODO: check that node id and location does not already exist
        self.nodes.append(Node(id, coord))

        # raise exception if duplicate added

    def add_element(self, element):
        """Adds an element to the fea class.

        Adds a FiniteElement object to the fea object. Checks for an exisiting
        id.

        Args:
            element:    FiniteElement object to be added to the fea object.

        Returns:
            void

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check that element id does not already exist
        self.elements.append(element)

        # raise exception if duplicate added

    def add_freedom_case(self, id, items=None):
        """Adds a freedom case to the fea class.

        Adds a FreedomCase object to the fea object. Checks for an exisiting
        id.

        Args:
            id:             FreedomCase id.
            items:          List of items to initialise the FreedomCase with.

        Returns:
            fc:             FreedomCase object.

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check to see if id already exists and raise exception

        fc = Cases.FreedomCase(self, id, items)
        self.freedom_cases.append(fc)
        return fc

    def add_load_case(self, id, items=None):
        """Adds a load case to the fea class.

        Adds a LoadCase object to the fea object. Checks for an exisiting id.

        Args:
            id:             LoadCase id.
            items:          List of items to initialise the LoadCase with.

        Returns:
            lc:             LoadCase object.

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check to see if id already exists and raise exception

        lc = Cases.LoadCase(self, id, items)
        self.load_cases.append(lc)
        return lc

    def add_analysis_case(self, id, fc_id, lc_id):
        """Adds an analysis case to the fea class.

        Adds an AnalysisCase object to the fea object. Checks for an exisiting
        id.

        Args:
            id:             AnalysisCase id.
            fc_id:          FreedomCase id to add to the analysis case.
            lc_id:          LoadCase id to add to the analysis case.

        Returns:
            void

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check to see if id already exists and raise exception

        self.analysis_cases.append(Cases.AnalysisCase(self, id, fc_id, lc_id))

    def find_node(self, node_id):
        """Finds the reference to the node object given its node_id.
        """

        # return the Node object whose id matches the given node_id
        for node in self.nodes:
            if node.id == node_id:
                return node
        else:
            raise FEAInputError("Cannot find node_id: {}".format(node_id))

    def get_node_lims(self):
        """asldkjasld
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
        """
        """

        # return the AnalysisCase object whose id matches the given case_id
        for analysis_case in self.analysis_cases:
            if analysis_case.id == case_id:
                return analysis_case
        else:
            raise FEAInputError("Cannot find AnalysisCase id: {}".format(
                case_id))


class FiniteElement:
    """
    """

    def __init__(self, analysis, id, node_ids):
        """Inits the FiniteElement class with an analysis object, element id
        and corresponding node ids.

        Args:
            analysis: reference to the analysis object.
            id: an integer representing a unique element id.
            node_ids: a list containing the node ids defining the element.

        Raises:
            TypeError: TODO
            ValueError: Raised if a negative id is provided. TODO
        """

        self.analysis = analysis
        self.id = id
        self.nodes = []
        self.f_int = ResultList()

        # add references to the node objects
        self.get_nodes(node_ids)

    def get_nodes(self, node_ids):
        """
        """

        for node_id in node_ids:
            try:
                self.nodes.append(self.analysis.find_node(node_id))
            except FEAInputError as error:
                print("Error in FiniteElement id: {}. {}".format(
                    self.id, error))

    def get_node_coords(self):
        """
        """

        coord_list = []

        for node in self.nodes:
            coord_list.append(node.coords)

        return np.array(coord_list)

    def get_nodal_displacements(self, case_id):
        """
        """

        disp_list = []  # allocate list of displacements

        for node in self.nodes:
            disp_list.append(node.get_displacement(case_id))

        return np.array(disp_list)

    def get_buckling_results(self, case_id, buckling_mode):
        """
        """

        v_list = []  # allocate list of eigenvectors

        for node in self.nodes:
            (w, v) = node.get_buckling_results(case_id, buckling_mode)
            v_list.append(v)

        return (w, np.array(v_list))

    def get_frequency_results(self, case_id, frequency_mode):
        """
        """

        v_list = []  # allocate list of eigenvectors

        for node in self.nodes:
            (w, v) = node.get_frequency_results(case_id, frequency_mode)
            v_list.append(v)

        return (w, np.array(v_list))

    def get_dofs(self):
        """
        """

        dof_list = np.array([], dtype=np.int32)

        for node in self.nodes:
            dof_list = np.append(dof_list, node.dofs)

        return dof_list

    def get_stiff_matrix(self):
        """
        """

        pass

    def set_fint(self, case_id, f_int):
        """
        """

        # save internal force vector in global coordinate system
        self.f_int.set_result(FrameForceVector(case_id, f_int))
