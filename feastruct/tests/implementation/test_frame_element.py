import unittest
import numpy as np
from feastruct.fea.node import Node
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
from feastruct.fea.elements.frame import FrameElement
from feastruct.fea.elements.frame2d import FrameElement2D
from feastruct.fea.cases import AnalysisCase


class TestFrameElement(unittest.TestCase):
    """Tests the functionality of the FrameElement class"""

    def setUp(self):
        steel = Steel()
        section = Section()

        # create an FrameElement2D elemnetn
        self.x = 4
        self.y = 3
        self.node1 = Node([0, 0])
        self.node2 = Node([self.x, self.y])

        self.element_frame = FrameElement(
            nodes=[self.node1, self.node2], material=steel, section=section,
            efs=[True, True, False, False, False, True]
        )
        self.element_frame2D = FrameElement2D(
            nodes=[self.node1, self.node2], material=steel, section=section,
            efs=[True, True, False, False, False, False]
        )

    def test_get_geometric_properties(self):
        (node_coords, dx, l0, c) = self.element_frame.get_geometric_properties()

        coord_list = np.array([
            [0, 0, 0],
            [self.x, self.y, 0],
        ])

        c_check = np.array([self.x, self.y, 0]) / np.sqrt(self.x ** 2 + self.y ** 2)

        self.assertEqual(node_coords.tolist(), coord_list.tolist())
        self.assertEqual(np.array([self.x, self.y, 0]).tolist(), dx.tolist())
        self.assertEqual(l0, np.sqrt(self.x ** 2 + self.y ** 2))
        self.assertEqual(c.tolist(), c_check.tolist())

    def test_get_nodal_displacements(self):
        # add dummy analysis case
        case1 = AnalysisCase('fc', 'lc')
        case2 = AnalysisCase('fc', 'lc')

        # create dummy displacements
        dof_count = 0

        for node in self.element_frame2D.nodes:
            for dof in node.dofs:
                dof.save_displacement(dof_count, case1)
                dof_count += 1

        dof_nums = [0, 1, 5]
        nodal_disps = self.element_frame2D.get_nodal_displacements(case1)

        for (i, node) in enumerate(nodal_disps):
            for (j, disp) in enumerate(node):
                self.assertEqual(disp, 6 * i + dof_nums[j])

        # check exception
        with self.assertRaises(Exception):
            self.element_frame2D.get_nodal_displacements(case2)


if __name__ == "__main__":
    unittest.main()
