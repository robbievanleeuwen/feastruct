import feastruct.fea.bcs as BCs
from feastruct.fea.exceptions import FEAInputError


class Case:
    """Parent class for cases.

    Provides an init method for the creation of different types of cases and
    also a method to add an item to the case list.

    :cvar analysis: Analysis object in which the case is used
    :vartype analysis: :class:`feastruct.fea.fea.fea`
    :cvar int id: A unique id that identifies the particular type of case
    :cvar items: A list containing entries in the case
    :vartype items: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
    """

    def __init__(self, analysis, id, items):
        """inits the Case class.

        :param analysis: Analysis object
        :type analysis: :class:`feastruct.fea.fea.fea`
        :param int id:  Unique id of the case
        :param items: A list of items to initialise the case with
        :type items: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
        """

        self.analysis = analysis
        self.id = id

        if items is None:
            self.items = []
        else:
            self.items = items

    def add_item(self, item):
        """Appends an 'item' to the list of entries in the case.

        :param item: Entry to add to the current case
        :type item: :class:`feastruct.fea.bcs.BoundaryCondition`
        """

        self.items.append(item)


class FreedomCase(Case):
    """Class for storing a set dirichlet boundary conditions.

    A freedom case contains a set of dirichlet boundary conditions that can be
    used in an analysis case. Methods are provided to add boundary conditions
    to the FreedomCase object.

    :cvar analysis: Analysis object in which the case is used
    :vartype analysis: :class:`feastruct.fea.fea.fea`
    :cvar int id: A unique id that identifies the particular type of case
    :cvar items: A list containing entries in the case
    :vartype items: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
    """

    def __init__(self, analysis, id, items):
        """inits the FreedomCase class.

        :param analysis: Analysis object
        :type analysis: :class:`feastruct.fea.fea.fea`
        :param int id:  Unique id of the case
        :param items: A list of items to initialise the case with
        :type items: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
        """

        super().__init__(analysis, id, items)

    def add_nodal_support(self, node_id, val, dir):
        """Adds a nodal dirichlet boundary condition to the current freedom
        case and to the node defined by the unique id 'node_id'.

        :param int node_id: Unique id of the node at which the boundary
            condition is applied
        :param float val: The value of the boundary condition
        :param int dir: The direction in which the boundary condition acts
        """

        # TODO: check that the support does not already exist
        # raise exception if duplicate added

        # add an entry to the freedom case items list
        self.add_item(BCs.NodalSupport(self.analysis, node_id, val, dir))


class LoadCase(Case):
    """Class for storing a set neumann boundary conditions.

    A load case contains a set of neumann boundary conditions that can be
    used in an analysis case. Methods are provided to add loads to the
    LoadCase object.

    :cvar analysis: Analysis object in which the case is used
    :vartype analysis: :class:`feastruct.fea.fea.fea`
    :cvar int id: A unique id that identifies the particular type of case
    :cvar items: A list containing entries in the case
    :vartype items: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
    """

    def __init__(self, analysis, id, items):
        """inits the LoadCase class.

        :param analysis: Analysis object
        :type analysis: :class:`feastruct.fea.fea.fea`
        :param int id:  Unique id of the case
        :param items: A list of items to initialise the case with
        :type items: list[:class:`feastruct.fea.bcs.BoundaryCondition`]
        """

        super().__init__(analysis, id, items)

    def add_nodal_load(self, node_id, val, dir):
        """Adds a nodal neumann boundary condition to the current load
        case and to the node defined by the unique id 'node_id'.

        :param int node_id: Unique id of the node at which the nodal load is
            applied
        :param float val: The value of the nodal load
        :param int dir: The direction in which the nodal load acts
        """

        # add an entry to the load case item list
        self.add_item(BCs.NodalLoad(self.analysis, node_id, val, dir))


class AnalysisCase:
    """Class for storing a specific freedom and load case to be used in an
    analysis case.

    An analysis case contains a reference to a combination of a freedom case
    and a load case, which are both used to define a particular analysis case.

    :cvar analysis: Analysis object in which the case is used
    :vartype analysis: :class:`feastruct.fea.fea.fea`
    :cvar int id: A unique id that identifies the particular type of case
    :cvar freedom_case: The FreedomCase object used in this analysis case
    :vartype freedom_case: :class:`feastruct.fea.cases.FreedomCase`
    :cvar load_case: The LoadCase object used in this analysis case
    :vartype load_case: :class:`feastruct.fea.cases.LoadCase`
    """

    def __init__(self, analysis, id, fc_id, lc_id):
        """inits the AnalysisCase class.

        :param analysis: Analysis object
        :type analysis: :class:`feastruct.fea.fea.fea`
        :param int id:  Unique id of the case
        :param int fc_id:  Unique id of freedom case used in this analysis case
        :param int lc_id:  Unique id of load case used in this analysis case
        :raises FEAInputError: If the unique id corresponding to the freedom
            case or load case cannot be located in the analysis object
        """

        self.analysis = analysis
        self.id = id

        try:
            # find freedom case
            for fc in self.analysis.freedom_cases:
                if fc.id == fc_id:
                    self.freedom_case = fc
                    break
            else:
                raise FEAInputError("Cannot find FreedomCase id: {}".format(
                    fc_id))

            # find load case
            for lc in self.analysis.load_cases:
                if lc.id == lc_id:
                    self.load_case = lc
                    break
            else:
                raise FEAInputError("Cannot find LoadCase id: {}".format(
                    lc_id))

        except FEAInputError as error:
            print(error)
