from feastruct.fea.fea import FiniteElementAnalysis
from feastruct.post.post2d import PostProcessor2D
import feastruct.fea.elements.frame2d as frame2d
import feastruct.fea.elements.frame3d as frame3d


class FrameAnalysis(FiniteElementAnalysis):
    """Parent class for a frame analysis.

    Includes a method for analysis initiation and a method for adding a frame element to the
    analysis.

    :cvar nodes: Nodes used in the finite element analysis
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis
    :vartype elements: list[:class:`~feastruct.fea.frame.FrameElement`]
    :cvar nfa: Node freedom arrangement
    :vartype nfa: list[bool]
    """

    def __init__(self, nodes, elements, nfa):
        """Inits the FrameAnalysis class.

        :param nodes: List of nodes with which to initialise the class
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class
        :type elements: list[:class:`~feastruct.fea.frame.FrameElement`]
        :param nfa: Node freedom arrangement
        :type nfa: list[bool]
        """

        # initialise parent FiniteElementAnalysis class
        super().__init__(nodes=nodes, elements=elements, nfa=nfa)

    def create_element(self, el_type, nodes, material, section):
        """Creates and returns a frame element and adds it to the
        :class:`~feastruct.fea.frame.Frame` object. Refer to 'xxx' for the possible element types.

        :param string el_type: String characterising the type of frame element
        :param nodes: List of node objects defining the element
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param material: Material object for the element
        :type material: :class:`~feastruct.pre.material.Material`
        :param section: Section object for the element
        :type section: :class:`~feastruct.pre.section.Section`
        :returns: FrameElement object
        :rtype: :class:`~feastruct.fea.frame.FrameElement`
        """

        if el_type == 'Bar2-2D':
            element = frame2d.Bar2D_2N(nodes=nodes, material=material, section=section)
        elif el_type == 'EB2-2D':
            element = frame2d.EulerBernoulli2D_2N(nodes=nodes, material=material, section=section)
        elif el_type == 'Bar2-3D':
            element = frame3d.Bar3D_2N(nodes=nodes, material=material, section=section)

        return(FiniteElementAnalysis.create_element(self, element=element))


class FrameAnalysis2D(FrameAnalysis):
    """Parent class for a 2D frame analysis.

    Includes a method for analysis initiation and a method for adding 2D frame elements to the
    analysis.

    Active degrees of freedom are *[x, y, rz]* and are referred to as *[0, 1, 5]*.

    :cvar nodes: Nodes used in the finite element analysis
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis
    :vartype elements: list[:class:`~feastruct.fea.frame.FrameElement`]
    :cvar nfa: Node freedom arrangement
    :vartype nfa: list[bool]
    :cvar post: Post-processor object
    :vartype post: :class:`feastruct.post.post2d.PostProcessor`
    """

    def __init__(self, nodes=None, elements=None):
        """Inits the FrameAnalysis2D class.

        :param nodes: List of nodes with which to initialise the class
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class
        :type elements: list[:class:`~feastruct.fea.frame.FrameElement`]
        """

        # set the node freedom arrangement for a 2D frame analysis
        nfa = [True, True, False, False, False, True]

        # initialise parent FrameAnalysis class
        super().__init__(nodes=nodes, elements=elements, nfa=nfa)

        self.post = PostProcessor2D(self)


class FrameAnalysis3D(FrameAnalysis):
    """Parent class for a 3D frame analysis.

    Includes a method for analysis initiation and a method for adding 3D frame elements to the
    analysis.

    Active degrees of freedom are *[x, y, z, rx, ry, rz]* and are referred to as
    *[0, 1, 2, 3, 4, 5]*.

    :cvar nodes: Nodes used in the finite element analysis
    :vartype nodes: list[:class:`~feastruct.fea.node.Node`]
    :cvar elements: Elements used in the finite element analysis
    :vartype elements: list[:class:`~feastruct.fea.frame.FrameElement`]
    :cvar nfa: Node freedom arrangement
    :vartype nfa: list[bool]
    """

    def __init__(self, nodes=None, elements=None):
        """Inits the FrameAnalysis3D class.

        :param nodes: List of nodes with which to initialise the class
        :type nodes: list[:class:`~feastruct.fea.node.Node`]
        :param elements: List of elements with which to initialise the class
        :type elements: list[:class:`~feastruct.fea.frame.FrameElement`]
        """

        # set the node freedom arrangement for a 3D frame analysis
        nfa = [True, True, True, True, True, True]

        # initialise parent FrameAnalysis class
        super().__init__(nodes=nodes, elements=elements, nfa=nfa)
