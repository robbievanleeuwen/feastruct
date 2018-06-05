import feastruct.fea.bcs as BCs
from feastruct.fea.exceptions import FEAInputError


class Case:
    """Parent class for cases.

    Provides an init method for the creation of different types of cases and
    also a method to add an item to the case list.

    Attributes:
        analysis:   An fea object in which the case is used [fea].
        id:         A unique id that identifies the particular type of
                    case [int].
        items:      A list containing entries in the case, e.g. nodal supports
                    or nodal loads.
    """

    def __init__(self, analysis, id, items):
        """inits the Case class.

        Args:
            analysis:   fea analysis object [fea]
            id:         Unique id of the case [int].
            items:      A list of items to initialise the case with.

        Returns:
            void
        """

        self.analysis = analysis
        self.id = id

        if items is None:
            self.items = []
        else:
            self.items = items

    def add_item(self, item):
        """Appends an 'item' to the list of entries in the case.

        Args:
            item:   Entry to add to the current case.

        Returns:
            void
        """

        self.items.append(item)


class FreedomCase(Case):
    """Class for storing a set dirichlet boundary conditions.

    A freedom case contains a set of dirichlet boundary conditions that can be
    used in an analysis case. Methods are provided to add boundary conditions
    to the FreedomCase object.

    Attributes:
        analysis:   An fea object in which the case is used [fea].
        id:         A unique id that identifies the particular type of
                    case [int].
        items:      A list containing entries in the case, e.g. nodal supports
                    or nodal loads.
    """

    def __init__(self, analysis, id, items):
        """inits the FreedomCase class.

        Args:
            analysis:   fea analysis object [fea]
            id:         Unique id of the case [int].
            items:      A list of items to initialise the case with.

        Returns:
            void
        """

        super().__init__(analysis, id, items)

    def add_nodal_support(self, node_id, val, dir):
        """Adds a nodal dirichlet boundary condition to the current freedom
        case and to the node defined by the unique id 'node_id'.

        Args:
            node_id:    Unique id of the node at which the boundary condition
                        is applied [int].
            val:        The value of the boundary condition [float].
            dir:        The direction in which the boundary condition
                        acts [int].

        Returns:
            void
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

    Attributes:
        analysis:   An fea object in which the case is used [fea].
        id:         A unique id that identifies the particular type of
                    case [int].
        items:      A list containing entries in the case, e.g. nodal supports
                    or nodal loads.
    """

    def __init__(self, analysis, id, items):
        """inits the LoadCase class.

        Args:
            analysis:   fea analysis object [fea]
            id:         Unique id of the case [int].
            items:      A list of items to initialise the case with.

        Returns:
            void
        """

        super().__init__(analysis, id, items)

    def add_nodal_load(self, node_id, val, dir):
        """Adds a nodal neumann boundary condition to the current load
        case and to the node defined by the unique id 'node_id'.

        Args:
            node_id:    Unique id of the node at which the boundary condition
                        is applied [int].
            val:        The value of the boundary condition [float].
            dir:        The direction in which the boundary condition
                        acts [int].

        Returns:
            void
        """

        # add an entry to the load case item list
        self.add_item(BCs.NodalLoad(self.analysis, node_id, val, dir))


class AnalysisCase:
    """Class for storing a specific freedom and load case to be used in an
    analysis case.

    An analysis case contains a reference to a combination of a freedom case
    and a load case, which are both used to define a particular analysis case.

    Attributes:
        analysis:       An fea object in which the case is used [fea].
        id:             A unique id that identifies the particular type of
                        case [int].
        freedom_case:   The FreedomCase object used in this analysis case
                        [FreedomCase].
        load_cae:       The LoadCase object used in this analysis case
                        [LoadCase].
    """

    def __init__(self, analysis, id, fc_id, lc_id):
        """inits the AnalysisCase class.

        Args:
            analysis:   fea analysis object [fea]
            id:         Unique id of the case [int].
            fc_id:      Unique id of freedom case used in this analysis
                        case [int].
            lc_id:      Unique id of load case used in this analysis
                        case [int].

        Returns:
            void

        Raises:
            FEAInputError:  If the unique id corresponding to the freedom case
                            or load case cannot be located in the analysis
                            object.
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
