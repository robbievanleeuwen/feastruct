import unittest
import numpy as np
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame import FrameAnalysis3D
from feastruct.solvers.linstatic import LinearStatic


class TestChapter5(unittest.TestCase):
    """Tests examples and problems from Chapter 5 of the mastan2 textbook."""

    def setUp(self):
        pass

    def test_example5_4(self):
        # create materials
        steel = Steel()

        # create 3d frame analysis object
        analysis = FrameAnalysis3D()

        section_ab = Section(area=20e3)
        section_ac = Section(area=30e3)
        section_ad = Section(area=40e3)
        section_ae = Section(area=30e3)

        # create nodes
        node_a = analysis.create_node(coords=[2000, 4000, 8000])
        node_b = analysis.create_node(coords=[0, 0, 0])
        node_c = analysis.create_node(coords=[8000, 0, 0])
        node_d = analysis.create_node(coords=[8000, 6000, 0])
        node_e = analysis.create_node(coords=[0, 6000, 0])

        # create truss elements
        element_ab = analysis.create_element(
            el_type='Bar2-3D', nodes=[node_a, node_b], material=steel, section=section_ab
        )
        element_ac = analysis.create_element(
            el_type='Bar2-3D', nodes=[node_a, node_c], material=steel, section=section_ac
        )
        element_ad = analysis.create_element(
            el_type='Bar2-3D', nodes=[node_a, node_d], material=steel, section=section_ad
        )
        element_ae = analysis.create_element(
            el_type='Bar2-3D', nodes=[node_a, node_e], material=steel, section=section_ae
        )

        # add supports
        freedom_case = cases.FreedomCase()
        sup_b = [0, 0, 0]
        sup_c = [0, 0, 0]
        sup_d = [0, 0, 0]
        sup_e = [0, 0, 0]
        sup_b[0] = freedom_case.add_nodal_support(node=node_b, val=0, dof=0)
        sup_b[1] = freedom_case.add_nodal_support(node=node_b, val=0, dof=1)
        sup_b[2] = freedom_case.add_nodal_support(node=node_b, val=0, dof=2)
        sup_c[0] = freedom_case.add_nodal_support(node=node_c, val=0, dof=0)
        sup_c[1] = freedom_case.add_nodal_support(node=node_c, val=0, dof=1)
        sup_c[2] = freedom_case.add_nodal_support(node=node_c, val=0, dof=2)
        sup_d[0] = freedom_case.add_nodal_support(node=node_d, val=0, dof=0)
        sup_d[1] = freedom_case.add_nodal_support(node=node_d, val=0, dof=1)
        sup_d[2] = freedom_case.add_nodal_support(node=node_d, val=0, dof=2)
        sup_e[0] = freedom_case.add_nodal_support(node=node_e, val=0, dof=0)
        sup_e[1] = freedom_case.add_nodal_support(node=node_e, val=0, dof=1)
        sup_e[2] = freedom_case.add_nodal_support(node=node_e, val=0, dof=2)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=node_a, val=200e3, dof=0)
        load_case.add_nodal_load(node=node_a, val=600e3, dof=1)
        load_case.add_nodal_load(node=node_a, val=-800e3, dof=2)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        dofs_a = node_a.get_dofs(node_a.nfs)
        u_a = dofs_a[0].get_displacement(analysis_case)
        v_a = dofs_a[1].get_displacement(analysis_case)
        w_a = dofs_a[2].get_displacement(analysis_case)

        self.assertEqual(np.around(u_a, 3), 0.178)
        self.assertEqual(np.around(v_a, 3), 2.722)
        self.assertEqual(np.around(w_a, 3), -0.487)

        # check reactions
        rb_x = sup_b[0].get_reaction(analysis_case=analysis_case)
        rb_y = sup_b[1].get_reaction(analysis_case=analysis_case)
        rb_z = sup_b[2].get_reaction(analysis_case=analysis_case)
        rc_x = sup_c[0].get_reaction(analysis_case=analysis_case)
        rc_y = sup_c[1].get_reaction(analysis_case=analysis_case)
        rc_z = sup_c[2].get_reaction(analysis_case=analysis_case)
        rd_x = sup_d[0].get_reaction(analysis_case=analysis_case)
        rd_y = sup_d[1].get_reaction(analysis_case=analysis_case)
        rd_z = sup_d[2].get_reaction(analysis_case=analysis_case)
        re_x = sup_e[0].get_reaction(analysis_case=analysis_case)
        re_y = sup_e[1].get_reaction(analysis_case=analysis_case)
        re_z = sup_e[2].get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(rb_x/1e3, 1), -76.4)
        self.assertEqual(np.around(rb_y/1e3, 1), -152.8)
        self.assertEqual(np.around(rb_z/1e3, 1), -305.6)
        self.assertEqual(np.around(rc_x/1e3, 1), 170.8)
        self.assertEqual(np.around(rc_y/1e3, 1), -113.9)
        self.assertEqual(np.around(rc_z/1e3, 1), -227.8)
        self.assertEqual(np.around(rd_x/1e3, 1), -470.8)
        self.assertEqual(np.around(rd_y/1e3, 1), -156.9)
        self.assertEqual(np.around(rd_z/1e3, 1), 627.8)
        self.assertEqual(np.around(re_x/1e3, 1), 176.4)
        self.assertEqual(np.around(re_y/1e3, 1), -176.4)
        self.assertEqual(np.around(re_z/1e3, 1), 705.6)


if __name__ == "__main__":
    unittest.main()
