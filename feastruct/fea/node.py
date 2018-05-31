from fea.exceptions import FEAInputError


class Node:
    """
    """

    def __init__(self, id, coord):
        """Inits the Node class with an id and coordinate.

        Args:
            id: an integer representing a unique node id.
            coord: a list consisting of the x, y coordinates of
            the node.

        Raises:
            TypeError: TODO
            ValueError: Raised if a negative id is provided. TODO
        """

        self.id = id
        self.coord = coord
        self.dofs = []
        self.u = []
        self.eigenvector = []
        self.fixity = [0, 0, 0]  # for post processing only
        # TODO: check value types

    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]

    @property
    def coords(self):
        return [self.x, self.y]

    def get_displacement(self, case_id):
        """
        """

        # get dictionary displacement entry for given case_id
        disp = next(d for d in self.u if d["case_id"] == case_id)

        if disp is not None:
            return disp["u"]
        else:
            raise FEAInputError("""Cannot find an analysis result for
            case_id: {} at node_id: {}""".format(case_id, self.id))

    def get_eigenvector(self, case_id, buckling_mode):
        """
        """

        # get dictionary eigenvector entry for given case_id
        v = next(d for d in self.eigenvector if d["case_id"] == case_id)

        if v is not None:
            return (v["v"][:, buckling_mode-1], v["w"][buckling_mode-1])
        else:
            raise FEAInputError("""Cannot find an eigenvector result for
            case_id: {} at node_id: {}""".format(case_id, self.id))
