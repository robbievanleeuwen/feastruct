import sys
import numpy as np

from fea.node import Node


class fea:
    """Summary of class here.

    Longer class information....
    Longer class information....

    Attributes:
        likes_spam: A boolean indicating if we like SPAM or not.
        eggs: An integer count of the eggs we have laid.
    """

    def __init__(self, analysis_type, nodes, elements):
        """Inits the fea class with..."""

        self.analysis_type = analysis_type
        self.dofs = 0
        self.nodes = nodes
        self.elements = elements
        self.supports = []
        self.nodal_loads = []

    def add_node(self, id, coord):
        """Adds a node to the fea class."""

        # TODO: check that node id and location does not already exist
        self.nodes.append(Node(id, coord))

        # raise exception if duplicate added

    def add_element(self, element):
        """Adds an element to the fea class."""

        # TODO: check that element id and nodes do not already exist
        self.elements.append(element)

        # raise exception if duplicate added

    def add_support(self, node_id, val, dir):
        """
        """

        # TODO: check that the support does not already exist
        self.supports.append({"node": node_id, "val": val, "dir": dir})

        # raise exception if duplicate added

    def add_nodal_load(self, node_id, val, dir):
        """
        """

        # TODO: check that the support does not already exist
        self.nodal_loads.append({"node": node_id, "val": val, "dir": dir})

        # raise exception if duplicate added

    def find_node(self, node_id):
        """Finds the reference to the node object given its node_id.
        """

        # return the node object whose id matches the given node_id
        for node in self.nodes:
            if node.id == node_id:
                return node
        else:
            raise IndexError("Cannot find node_id: {}".format(node_id))


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

        # add references to the node objects
        self.get_nodes(node_ids)

    def get_nodes(self, node_ids):
        """
        """

        for node_id in node_ids:
            try:
                self.nodes.append(self.analysis.find_node(node_id))
            except IndexError as error:
                print(
                    "Error in FiniteElement id: {}. {}".format(self.id, error))
                sys.exit(1)

    def get_coords(self):
        """
        """

        coord_list = []

        for node in self.nodes:
            coord_list.append(node.coords)

        return coord_list

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
