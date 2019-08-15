class ScalarResult:
    """Class for storing a scalar result for a specific analysis case.

    :cvar float result: Scalar result
    :cvar analysis_case: Analysis case relating to the result
    :vartype analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    """

    def __init__(self, result, analysis_case):
        """Inits the ScalarResult class.

        :cvar float result: Scalar result
        :param analysis_case: Analysis case relating to the result
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        self.result = result
        self.analysis_case = analysis_case
