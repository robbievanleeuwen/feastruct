import feastruct.post.results as results


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
        self.fixity = [0, 0, 0]  # for post processing only
        self.u = results.ResultList()
        self.buckling_v = results.ResultList()
        self.frequency_v = results.ResultList()
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

        return self.u.get_result(case_id).u

    def set_displacements(self, case_id, u):
        """
        """

        self.u.set_result(results.Displacement(case_id, u))

    def get_buckling_results(self, case_id, buckling_mode):
        """
        """

        result = self.buckling_v.get_result(case_id, mode=buckling_mode)

        return (result.w, result.v)

    def set_buckling_results(self, case_id, buckling_mode, w, v):
        """
        """

        self.buckling_v.set_result(
            results.EigenResult(case_id, buckling_mode, w, v),
            mode=buckling_mode)

    def get_frequency_results(self, case_id, frequency_mode):
        """
        """

        result = self.frequency_v.get_result(case_id, mode=frequency_mode)

        return (result.w, result.v)

    def set_frequency_results(self, case_id, frequency_mode, w, v):
        """
        """

        self.frequency_v.set_result(
            results.EigenResult(case_id, frequency_mode, w, v),
            mode=frequency_mode)
