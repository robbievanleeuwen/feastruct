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
        self.fixity = [0, 0, 0]  # for post processing only

    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]

    @property
    def coords(self):
        return [self.x, self.y]
