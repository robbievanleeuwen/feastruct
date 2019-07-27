import unittest
import numpy as np
from feastruct.fea.node import Node
from feastruct.fea.cases import AnalysisCase


class TestNode(unittest.TestCase):
    """Tests the functionality of the Node class"""

    def setUp(self):
        self.x = 5
        self.y = 10
        self.z = 0
        self.node = Node([self.x, self.y, self.z])

    def test_coordinates(self):
        self.assertTrue(np.isclose(self.node.x, self.x))
        self.assertTrue(np.isclose(self.node.y, self.y))
        self.assertTrue(np.isclose(self.node.z, self.z))

    def test_dof_size(self):
        self.assertTrue(len(self.node.dofs) == 6)

    def test_dof_parent(self):
        for dof in self.node.dofs:
            self.assertTrue(dof.node is self.node)

    def test_node_dof_nums(self):
        for (i, dof) in enumerate(self.node.dofs):
            self.assertTrue(dof.node_dof_num is i)

    def test_node_get_dofs(self):
        # case 1 - 3D dofs
        node_dof_nums = [0, 1, 2, 3, 4, 5]

        dof_list = self.node.get_dofs(node_dof_nums)

        for (i, dof) in enumerate(dof_list):
            self.assertTrue(dof.node_dof_num is node_dof_nums[i])

        # case 2 - 2D dofs
        node_dof_nums = [0, 1, 3]

        dof_list = self.node.get_dofs(node_dof_nums)

        for (i, dof) in enumerate(dof_list):
            self.assertTrue(dof.node_dof_num is node_dof_nums[i])

        # case 3 - single DoF
        node_dof_nums = [4]

        dof_list = self.node.get_dofs(node_dof_nums)

        for (i, dof) in enumerate(dof_list):
            self.assertTrue(dof.node_dof_num is node_dof_nums[i])

    def test_dof_displacements(self):
        # add dummy analysis case
        case1 = AnalysisCase('fc', 'lc')
        case2 = AnalysisCase('fc', 'lc')

        # create dummy displacements
        for (i, dof) in enumerate(self.node.dofs):
            dof.save_displacement(i, case1)

        # get displacements
        for (i, dof) in enumerate(self.node.dofs):
            self.assertEqual(dof.get_displacement(case1), i)

        # check exception
        with self.assertRaises(Exception):
            dof.get_displacement(case2)


if __name__ == "__main__":
    unittest.main()
