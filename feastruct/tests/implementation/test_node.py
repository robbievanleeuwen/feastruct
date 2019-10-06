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
        freedom_signature = [True, True, True, True, True, True]
        indices = [i for i, x in enumerate(freedom_signature) if x]

        dof_list = self.node.get_dofs(freedom_signature)

        for (i, dof) in enumerate(dof_list):
            self.assertTrue(dof.node_dof_num == indices[i])

        # case 2 - 2D dofs
        freedom_signature = [True, True, False, False, False, True]
        indices = [i for i, x in enumerate(freedom_signature) if x]

        dof_list = self.node.get_dofs(freedom_signature)

        for (i, dof) in enumerate(dof_list):
            self.assertTrue(dof.node_dof_num == indices[i])

        # case 3 - single DoF
        freedom_signature = [True, False, False, False, False, False]
        indices = [i for i, x in enumerate(freedom_signature) if x]

        dof_list = self.node.get_dofs(freedom_signature)

        for (i, dof) in enumerate(dof_list):
            self.assertTrue(dof.node_dof_num == indices[i])

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

    def test_node_displacements(self):
        # add dummy analysis case
        case1 = AnalysisCase('fc', 'lc')

        # create dummy displacements
        for (i, dof) in enumerate(self.node.dofs):
            # don't save dof index 4
            if i == 4:
                pass
            else:
                dof.save_displacement(i, case1)

        # get displacements
        disp_vector = self.node.get_displacements(case1)

        for (i, disp) in enumerate(disp_vector):
            if i == 4:
                self.assertEqual(disp, None)
            else:
                self.assertEqual(disp, i)


if __name__ == "__main__":
    unittest.main()
