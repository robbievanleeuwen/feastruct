from feastruct.fea.exceptions import FEAInputError


class ResultList:
    """
    """

    def __init__(self):
        """"""

        self.results = []

    def get_result(self, case_id, mode=None):
        """"""
        # find case_id in results list
        for result in self.results:
            if result.case_id == case_id:
                if mode is None:
                    return result
                else:
                    if result.mode == mode:
                        return result

        # if case_id cannot be found
        raise FEAInputError("Cannot find a result for case_id: {}".format(
            case_id))

    def set_result(self, new_result, mode=None):
        """"""

        # check to see if case_id already exists
        for result in self.results:
            if result.case_id == new_result.case_id:
                if mode is None:
                    # replace result if case_id already exists
                    result = new_result
                    return
                else:
                    if result.mode == mode:
                        # replace result if case_id and mode already exists
                        result = new_result
                        return
                    else:
                        # append to result list if this is a new case_id
                        self.results.append(new_result)

        # append to result list if this is a new case_id
        self.results.append(new_result)


class ResultItem:
    """
    """

    def __init__(self, case_id):
        self.case_id = case_id


class Displacement(ResultItem):
    """
    """

    def __init__(self, case_id, u):
        """"""

        super().__init__(case_id)
        self.u = u


class EigenResult(ResultItem):
    """
    """

    def __init__(self, case_id, mode, w, v):
        """"""

        super().__init__(case_id)
        self.mode = mode
        self.w = w
        self.v = v


class Force(ResultItem):
    """
    """

    def __init__(self, case_id, f):
        """"""

        super().__init__(case_id)
        self.f = f


class FrameForceVector((ResultItem)):
    """
    """

    def __init__(self, case_id, f):
        """"""

        super().__init__(case_id)
        self.N1 = f[0]
        self.V1 = f[1]
        self.M1 = f[2]
        self.N2 = f[3]
        self.V2 = f[4]
        self.M2 = f[5]
