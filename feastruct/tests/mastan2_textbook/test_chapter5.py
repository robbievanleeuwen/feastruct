import unittest
import numpy as np
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D, FrameAnalysis3D
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
        disp_a = node_a.get_displacements(analysis_case)
        u_a = disp_a[0]
        v_a = disp_a[1]
        w_a = disp_a[2]

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

    def test_example5_6(self):
        # create materials
        steel = Steel()

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        section_ab = Section(area=6e3, ixx=200e6)
        section_bc = Section(area=4e3, ixx=50e6)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[8000])
        node_p = analysis.create_node(coords=[10000])
        node_c = analysis.create_node(coords=[13000])

        # create beam elements
        element_ab = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=steel, section=section_ab
        )
        element_bp = analysis.create_element(
            el_type='EB2-2D', nodes=[node_b, node_p], material=steel, section=section_bc
        )
        element_pc = analysis.create_element(
            el_type='EB2-2D', nodes=[node_p, node_c], material=steel, section=section_bc
        )

        # add supports
        freedom_case = cases.FreedomCase()
        sup_a = freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        sup_b = freedom_case.add_nodal_support(node=node_b, val=0, dof=1)
        sup_c = [0, 0, 0]
        sup_c[0] = freedom_case.add_nodal_support(node=node_c, val=0, dof=0)
        sup_c[1] = freedom_case.add_nodal_support(node=node_c, val=0, dof=1)
        sup_c[2] = freedom_case.add_nodal_support(node=node_c, val=0, dof=5)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_element_load(element_ab.generate_udl(q=-2))
        load_case.add_nodal_load(node=node_p, val=-20e3, dof=1)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_a = node_a.get_displacements(analysis_case)
        disp_b = node_b.get_displacements(analysis_case)
        r_a = disp_a[5]
        r_b = disp_b[5]

        self.assertEqual(np.around(r_a, 7), -5.681e-4)
        self.assertEqual(np.around(r_b, 7), 0.696e-4)

        # check reactions
        ra_y = sup_a.get_reaction(analysis_case=analysis_case)
        rb_y = sup_b.get_reaction(analysis_case=analysis_case)
        rc_y = sup_c[1].get_reaction(analysis_case=analysis_case)
        rc_m = sup_c[2].get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(ra_y/1e3, 2), 6.13)
        self.assertEqual(np.around(rb_y/1e3, 2), 23.00)
        self.assertEqual(np.around(rc_y/1e3, 2), 6.87)
        self.assertEqual(np.around(rc_m/1e6, 2), -9.32)

    def test_example5_7(self):
        l_a = np.sqrt(8000 * 8000 - 3000 * 3000)

        # create materials
        steel = Steel()

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        section = Section(area=6e3, ixx=200e6)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[l_a, 3000])
        node_c = analysis.create_node(coords=[l_a + 8000, 3000])

        # create beam elements
        element_ab = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=steel, section=section
        )
        element_bc = analysis.create_element(
            el_type='EB2-2D', nodes=[node_b, node_c], material=steel, section=section
        )

        # add supports
        freedom_case = cases.FreedomCase()
        sup_a = [0, 0, 0]
        sup_a[0] = freedom_case.add_nodal_support(node=node_a, val=0, dof=0)
        sup_a[1] = freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        sup_a[2] = freedom_case.add_nodal_support(node=node_a, val=0, dof=5)
        sup_c = [0, 0, 0]
        sup_c[0] = freedom_case.add_nodal_support(node=node_c, val=0, dof=0)
        sup_c[1] = freedom_case.add_nodal_support(node=node_c, val=0, dof=1)
        sup_c[2] = freedom_case.add_nodal_support(node=node_c, val=0, dof=5)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=node_b, val=50e3 * 3000 / 8000, dof=0)
        load_case.add_nodal_load(node=node_b, val=-50e3 * l_a / 8000, dof=1)
        load_case.add_element_load(element_bc.generate_udl(q=-4))

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_b = node_b.get_displacements(analysis_case)
        u_b = disp_b[0]
        v_b = disp_b[1]
        r_b = disp_b[5]

        self.assertEqual(np.around(u_b, 4), 0.9950)
        self.assertEqual(np.around(v_b, 3), -4.982)
        self.assertEqual(np.around(r_b, 6), -0.000534)

        # check reactions
        ra_x = sup_a[0].get_reaction(analysis_case=analysis_case)
        ra_y = sup_a[1].get_reaction(analysis_case=analysis_case)
        ra_m = sup_a[2].get_reaction(analysis_case=analysis_case)
        rc_x = sup_c[0].get_reaction(analysis_case=analysis_case)
        rc_y = sup_c[1].get_reaction(analysis_case=analysis_case)
        rc_m = sup_c[2].get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(ra_x/1e3, 1), 130.5)
        self.assertEqual(np.around(ra_y/1e3, 1), 55.7)
        self.assertEqual(np.around(ra_m/1e6, 2), 13.37)
        self.assertEqual(np.around(rc_x/1e3, 1), -149.3)
        self.assertEqual(np.around(rc_y/1e3, 1), 22.7)
        self.assertEqual(np.around(rc_m/1e6, 2), -45.36)


if __name__ == "__main__":
    unittest.main()
