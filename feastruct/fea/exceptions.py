class FEAInputError(Exception):
    """Raise for errors related to pre-processing."""


class FEASolverError(Exception):
    """Raise for errors related to FEA solvers."""


class FEAPostError(Exception):
    """Raise for errors related to FEA post-processing."""
