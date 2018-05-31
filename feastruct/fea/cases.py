import sys
import fea.bcs as BCs
from fea.exceptions import FEAInputError


class Case:
    """
    """

    def __init__(self, analysis, id, items):
        """"""

        self.analysis = analysis
        self.id = id

        if items is None:
            self.items = []
        else:
            self.items = items

    def add_item(self, item):
        """"""

        self.items.append(item)


class FreedomCase(Case):
    """
    """

    def __init__(self, analysis, id, items):
        super().__init__(analysis, id, items)

    def add_nodal_support(self, node_id, val, dir):
        """
        """

        # TODO: check that the support does not already exist
        # raise exception if duplicate added

        # add a dictionary entry to the supports list
        self.add_item(BCs.NodalSupport(self.analysis, node_id, val, dir))


class LoadCase(Case):
    """
    """

    def __init__(self, analysis, id, items):
        super().__init__(analysis, id, items)

    def add_nodal_load(self, node_id, val, dir):
        """
        """

        # TODO: check that the load does not already exist
        # raise exception if duplicate added

        # add a dictionary entry to the supports list
        self.add_item(BCs.NodalLoad(self.analysis, node_id, val, dir))


class AnalysisCase:
    """
    """

    def __init__(self, analysis, id, fc_id, lc_id):
        """
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
            sys.exit(1)
