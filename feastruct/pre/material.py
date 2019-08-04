class Material:
    """Class for structural materials (assumed isotropic).

    Provides a way of storing material properties related to a specific material. The colour can be
    a multitude of different formats, refer to https://matplotlib.org/api/colors_api.html and
    https://matplotlib.org/examples/color/named_colors.html for more information.

    :cvar string name: Material name
    :cvar float elastic_modulus: Material modulus of elasticity
    :cvar float poissons_ratio: Material Poisson's ratio
    :cvar float shear_modulus: Material shear modulus, derived from the elastic modulus and
        Poisson's ratio assuming an isotropic material
    :cvar float rho: Material density
    :cvar colour: Material colour for rendering
    :vartype colour: :class:`matplotlib.colors`
    """

    def __init__(self, name, elastic_modulus, poissons_ratio, rho, colour='w'):
        """Inits the Material class.

        :param string name: Material name
        :param float elastic_modulus: Material modulus of elasticity
        :param float poissons_ratio: Material Poisson's ratio
        :param float rho: Material density
        :param colour: Material color for rendering
        :type colour: :class:`matplotlib.colors`
        """

        self.name = name
        self.elastic_modulus = elastic_modulus
        self.poissons_ratio = poissons_ratio
        self.shear_modulus = elastic_modulus / (2 * (1 + poissons_ratio))
        self.rho = rho
        self.colour = colour


class Steel(Material):
    """Class for structural steel to AS4100-1998.

    Material properties:

    * name: steel
    * elastic_modulus: 200000 N/mm2
    * poissons_ratio: 0.3
    * rho: 7.85e-9 T/mm^3
    * colour: lightgrey
    """

    def __init__(self):
        """Inits the Steel class."""

        super().__init__(name='steel', elastic_modulus=2e5, poissons_ratio=0.3, rho=7.85e-9,
                         colour='lightgrey')
