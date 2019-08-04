import unittest
import numpy as np
from feastruct.fea.node import Node
from feastruct.pre.material import Steel
from feastruct.fea.fea import FiniteElement
from feastruct.fea.cases import AnalysisCase


class TestFiniteElement(unittest.TestCase):
    """Tests the functionality of the FiniteElement class"""

    def setUp(self):
        steel = Steel()

        # create a rectangular element
        self.width = 5
        self.height = 3

        self.node1 = Node([0, 0])
        self.node2 = Node([self.width, 0])
        self.node3 = Node([self.width, self.height])
        self.node4 = Node([0, self.height])

        self.element = FiniteElement(
            nodes=[self.node1, self.node2, self.node3, self.node4],
            material=steel, efs=[True, True, False, False, False, False]
        )

    def test_get_node_coords(self):
        coord_list = np.array([
            [0, 0, 0],
            [self.width, 0, 0],
            [self.width, self.height, 0],
            [0, self.height, 0]
        ])

        self.assertEqual(coord_list.tolist(), self.element.get_node_coords().tolist())

    def test_get_dofs(self):
        indices = [i for i, x in enumerate(self.element.efs) if x]

        dof_list = self.element.get_dofs()

        for (i, node) in enumerate(dof_list):
            for (j, dof) in enumerate(node):
                self.assertTrue(dof is self.element.nodes[i].dofs[indices[j]])

    def test_get_nodal_displacements(self):
        # add dummy analysis case
        case1 = AnalysisCase('fc', 'lc')
        case2 = AnalysisCase('fc', 'lc')

        # create dummy displacements
        dof_count = 0

        for node in self.element.nodes:
            for dof in node.dofs:
                dof.save_displacement(dof_count, case1)
                dof_count += 1

        indices = [i for i, x in enumerate(self.element.efs) if x]
        nodal_disps = self.element.get_nodal_displacements(case1)

        for (i, node) in enumerate(nodal_disps):
            for (j, disp) in enumerate(node):
                self.assertEqual(disp, 6 * i + indices[j])

        # check exception
        with self.assertRaises(Exception):
            self.element.get_nodal_displacements(case2)

    def test_fint(self):
        # add dummy analysis case
        case1 = AnalysisCase('fc', 'lc')
        case2 = AnalysisCase('fc', 'lc')
        case3 = AnalysisCase('fc', 'lc')

        # create dummy force vector
        f1 = np.array(range(8))
        f2 = 3 * np.array(range(8))

        # save fint
        self.element.save_fint(f1, case1)
        self.element.save_fint(f2, case2)

        # get fint
        f_res1 = self.element.get_fint(case1)
        f_res2 = self.element.get_fint(case2)

        self.assertEqual(f1.tolist(), f_res1.tolist())
        self.assertEqual(f2.tolist(), f_res2.tolist())

        # check exception
        with self.assertRaises(Exception):
            self.element.get_fint(case3)


if __name__ == "__main__":
    unittest.main()
