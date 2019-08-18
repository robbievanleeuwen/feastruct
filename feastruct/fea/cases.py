import feastruct.fea.bcs as bcs


class AnalysisCase:
    """Class for storing a specific freedom and load case to be used in an analysis case.

    An analysis case contains a reference to a specific freedom and load case, which are both used
    to define a particular analysis case. Results are stored making reference to an AnalysisCase
    object.

    :cvar freedom_case: The FreedomCase object used in this analysis case
    :vartype freedom_case: :class:`~feastruct.fea.cases.FreedomCase`
    :cvar load_case: The LoadCase object used in this analysis case
    :vartype load_case: :class:`~feastruct.fea.cases.LoadCase`
    """

    def __init__(self, freedom_case, load_case):
        """Inits the AnalysisCase class.

        :param freedom_case: Freedom case used in this analysis case
        :param load_case: Load case used in this analysis case
        """

        self.freedom_case = freedom_case
        self.load_case = load_case


class Case:
    """Parent class for cases.

    Provides an init method for the creation of different types of cases and also a method to add
    an item to the case list.

    :cvar items: A list of BoundaryCondition items defining the case
    :vartype items: list[:class:`~feastruct.fea.bcs.BoundaryCondition`]
    """

    def __init__(self, items):
        """Inits the Case class.

        :param items: A list of BoundaryConditions to initialise the case
        :type items: list[:class:`~feastruct.fea.bcs.BoundaryCondition`]
        """

        if items is None:
            self.items = []
        else:
            self.items = items

    def add_item(self, item):
        """Appends an 'item' to the list of entries in the case.

        :param item: Entry to add to the current case
        :type item: :class:`~feastruct.fea.bcs.BoundaryCondition`
        """

        self.items.append(item)


class FreedomCase(Case):
    """Class for storing a set dirichlet boundary conditions.

    A freedom case contains a set of dirichlet boundary conditions that can be used in an analysis
    case. Methods are provided to add boundary conditions to the FreedomCase object.

    :cvar items: A list of BoundaryCondition items defining the freedom case
    :vartype items: list[:class:`~feastruct.fea.bcs.BoundaryCondition`]
    """

    def __init__(self, items=None):
        """Inits the FreedomCase class.

        :param items: A list of BoundaryConditions to initialise the case
        :type items: list[:class:`~feastruct.fea.bcs.BoundaryCondition`]
        """

        super().__init__(items)

    def add_nodal_support(self, node, val, dof):
        """Adds a nodal dirichlet boundary condition to the current freedom case.

        :param node: The node object at which the nodal support acts
        :type node: :class:`~feastruct.fea.node.Node`
        :param float val: The value of the nodal support - zero indicates a fixed support, whereas
            a non-zero value indicates a prescribed displacement
        :param int dof: The degree of freedom about which the boundary condition acts

        :returns: Nodal support object
        :rtype: :class:`feastruct.fea.bcs.NodalSupport`
        """

        # TODO: check that the support does not already exist
        # raise exception if duplicate added

        # add an entry to the freedom case items list
        new_support = bcs.NodalSupport(node, val, dof)
        self.add_item(new_support)

        return new_support

    def get_nodal_fixities(self, node):
        """Returns a list defining the nodal fixity at the node for the current freedom case.

        :param node: Node object
        :type node: :class:`feastruct.fea.node.Node`
        """

        fixity = [0, 0, 0]

        # loop through all supports
        for support in self.items:
            # if the support is the node in question and the support is fixed
            if support.val == 0 and support.node == node:
                if support.dof == 5:
                    dof = 2
                else:
                    dof = support.dof

                fixity[dof] = 1

        return fixity


class LoadCase(Case):
    """Class for storing a set neumann boundary conditions.

    A load case contains a set of neumann boundary conditions that can be used in an analysis case.
    Methods are provided to add loads to the LoadCase object.

    :cvar items: A list of BoundaryCondition items defining the load case
    :vartype items: list[:class:`~feastruct.fea.bcs.BoundaryCondition`]
    :cvar element_items: A list of ElementLoad items defining the load case
    :vartype items: list[:class:`~feastruct.fea.bcs.ElementLoad`]
    """

    def __init__(self, items=None):
        """Inits the LoadCase class.

        :param items: A list of BoundaryConditions to initialise the case
        :type items: list[:class:`~feastruct.fea.bcs.BoundaryCondition`]
        """

        super().__init__(items)
        self.element_items = []

    def add_nodal_load(self, node, val, dof):
        """Adds a nodal neumann boundary condition to the current load case.

        :param node: The node object at which the nodal load is applied
        :type node: :class:`~feastruct.fea.node.Node`
        :param float val: The value of the nodal load
        :param int dof: The degree of freedom about which the nodal load acts

        :returns: Nodal load object
        :rtype: :class:`feastruct.fea.bcs.NodalLoad`
        """

        # add an entry to the load case item list
        new_load = bcs.NodalLoad(node, val, dof)
        self.add_item(new_load)

        return new_load

    def add_element_load(self, element_load):
        """Adds an element load to the current load case.

        :param element_load: Element load object
        :type element_load: :class:`~feastruct.fea.bcs.ElementLoad`

        :returns: Element load object
        :rtype: :class:`~feastruct.fea.bcs.ElementLoad`
        """

        # add an entry to the load case element_items list
        self.element_items.append(element_load)

        return element_load
