class Section:
    """Class for storing cross-section properties.

    :cvar string name: Section name
    :cvar float area: Cross-sectional area
    :cvar float ixx: Second moment of area about the centroidal x-axis
    :cvar float iyy: Second moment of area about the centroidal y-axis
    :cvar float j: Torsion constant
    :cvar float A_sx: Shear area about the x-axis
    :cvar float A_sy: Shear area about the y-axis
    :cvar colour: Section colour for rendering
    :vartype colour: :class:`matplotlib.colors`
    """

    def __init__(self, name='', area=1, ixx=1, iyy=1, j=1, A_sx=1, A_sy=1, colour='k'):
        """Inits the Section class.

        :param string name: Section name
        :param float area: Cross-sectional area
        :param float ixx: Second moment of area about the centroidal x-axis
        :param float iyy: Second moment of area about the centroidal y-axis
        :param float j: Torsion constant
        :param float A_sx: Shear area about the x-axis
        :param float A_sy: Shear area about the y-axis
        :param colour: Section colour for rendering
        :type colour: :class:`matplotlib.colors`
        """

        self.name = name
        self.area = area
        self.ixx = ixx
        self.iyy = iyy
        self.j = j
        self.A_sx = A_sx
        self.A_sy = A_sy
        self.colour = colour
