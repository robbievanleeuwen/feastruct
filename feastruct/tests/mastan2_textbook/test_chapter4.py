import unittest
import numpy as np
from feastruct.pre.material import Material, Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic


class TestChapter4(unittest.TestCase):
    """Tests examples and problems from Chapter 4 of the mastan2 textbook."""

    def setUp(self):
        # create materials
        steel = Steel()

        # create 2d frame analysis object
        self.analysis = FrameAnalysis2D()

        section1 = Section(area=6e3, ixx=200e6)
        section2 = Section(area=4e3, ixx=50e6)

        # create nodes
        self.node_a = self.analysis.create_node(coords=[0])
        self.node_b = self.analysis.create_node(coords=[8000])
        self.node_c = self.analysis.create_node(coords=[13000])

        # create beam elements
        self.element_ab = self.analysis.create_element(
            el_type='EB2-2D', nodes=[self.node_a, self.node_b], material=steel, section=section1
        )
        self.element_bc = self.analysis.create_element(
            el_type='EB2-2D', nodes=[self.node_b, self.node_c], material=steel, section=section2
        )

    def test_example4_9(self):
        # add supports
        freedom_case = cases.FreedomCase()
        sup_a = [0, 0, 0]
        sup_a[0] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=0)
        sup_a[1] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=1)
        sup_a[2] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=5)
        sup_b = freedom_case.add_nodal_support(node=self.node_b, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=self.node_c, val=5e3/np.sqrt(2), dof=0)
        load_case.add_nodal_load(node=self.node_c, val=-5e3/np.sqrt(2), dof=1)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=self.analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_b = self.node_b.get_displacements(analysis_case)
        disp_c = self.node_c.get_displacements(analysis_case)

        u_b = disp_b[0]
        rz_b = disp_b[5]
        u_c = disp_c[0]
        v_c = disp_c[1]
        rz_c = disp_c[5]

        self.assertEqual(np.around(u_b, 3), 0.024)
        self.assertEqual(np.around(rz_b, 5), -0.00088)
        self.assertEqual(np.around(u_c, 3), 0.046)
        self.assertEqual(np.around(v_c, 2), -19.15)
        self.assertEqual(np.around(rz_c, 4), -0.0053)

        # check reactions
        ra_x = sup_a[0].get_reaction(analysis_case=analysis_case)
        ra_y = sup_a[1].get_reaction(analysis_case=analysis_case)
        ra_m = sup_a[2].get_reaction(analysis_case=analysis_case)
        rb_y = sup_b.get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(ra_x/1e3, 2), -3.54)
        self.assertEqual(np.around(ra_y/1e3, 2), -3.31)
        self.assertEqual(np.around(ra_m/1e6, 2), -8.84)
        self.assertEqual(np.around(rb_y/1e3, 2), 6.85)

        # check bending moment forces
        (_, m_ab) = self.element_ab.get_bmd(n=2, analysis_case=analysis_case)
        (_, m_bc) = self.element_bc.get_bmd(n=2, analysis_case=analysis_case)

        self.assertTrue(np.isclose(m_ab[1], m_bc[0]))
        self.assertEqual(np.around(m_ab[1]/1e6, 2), 17.68)

    def test_example4_10(self):
        # add supports
        freedom_case = cases.FreedomCase()
        sup_a = [0, 0, 0]
        sup_c = [0, 0, 0]
        sup_a[0] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=0)
        sup_a[1] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=1)
        sup_a[2] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=5)
        sup_b = freedom_case.add_nodal_support(node=self.node_b, val=0, dof=1)
        sup_c[0] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=0)
        sup_c[1] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=1)
        sup_c[2] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=5)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=self.node_b, val=50e6, dof=5)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=self.analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_b = self.node_b.get_displacements(analysis_case)
        rz_b = disp_b[5]

        self.assertEqual(np.around(rz_b, 5), 0.00179)

        # check reactions
        ra_y = sup_a[1].get_reaction(analysis_case=analysis_case)
        ra_m = sup_a[2].get_reaction(analysis_case=analysis_case)
        rb_y = sup_b.get_reaction(analysis_case=analysis_case)
        rc_y = sup_c[1].get_reaction(analysis_case=analysis_case)
        rc_m = sup_c[2].get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(ra_y/1e3, 2), 6.70)
        self.assertEqual(np.around(ra_m/1e6, 2), 17.86)
        self.assertEqual(np.around(rb_y/1e3, 2), -2.41)
        self.assertEqual(np.around(rc_y/1e3, 2), -4.29)
        self.assertEqual(np.around(rc_m/1e6, 2), 7.14)

        # check bending moment forces
        (_, m_ab) = self.element_ab.get_bmd(n=2, analysis_case=analysis_case)
        (_, m_bc) = self.element_bc.get_bmd(n=2, analysis_case=analysis_case)

        self.assertEqual(np.around(m_ab[1]/1e6, 2), -35.71)
        self.assertEqual(np.around(m_bc[0]/1e6, 2), 14.29)

    def test_example4_11(self):
        # add supports
        freedom_case = cases.FreedomCase()
        sup_a = [0, 0, 0]
        sup_c = [0, 0, 0]
        sup_a[0] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=0)
        sup_a[1] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=1)
        sup_a[2] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=5)
        sup_b = freedom_case.add_nodal_support(node=self.node_b, val=-20, dof=1)
        sup_c[0] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=0)
        sup_c[1] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=1)
        sup_c[2] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=5)

        # add loads
        load_case = cases.LoadCase()

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=self.analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_b = self.node_b.get_displacements(analysis_case)
        rz_b = disp_b[5]

        self.assertEqual(np.around(rz_b, 6), -0.000964)

        # check reactions
        ra_y = sup_a[1].get_reaction(analysis_case=analysis_case)
        ra_m = sup_a[2].get_reaction(analysis_case=analysis_case)
        rb_y = sup_b.get_reaction(analysis_case=analysis_case)
        rc_y = sup_c[1].get_reaction(analysis_case=analysis_case)
        rc_m = sup_c[2].get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(ra_y/1e3, 2), 15.13)
        self.assertEqual(np.around(ra_m/1e6, 2), 65.36)
        self.assertEqual(np.around(rb_y/1e3, 2), -36.65)
        self.assertEqual(np.around(rc_y/1e3, 2), 21.51)
        self.assertEqual(np.around(rc_m/1e6, 2), -51.86)

        # check bending moment forces
        (_, m_ab) = self.element_ab.get_bmd(n=2, analysis_case=analysis_case)
        (_, m_bc) = self.element_bc.get_bmd(n=2, analysis_case=analysis_case)

        self.assertTrue(np.isclose(m_ab[1], m_bc[0]))
        self.assertEqual(np.around(m_ab[1]/1e6, 2), -55.71)

    def test_example4_13(self):
        # move node c
        self.node_c.coords = [8000, -5000, 0]

        # add supports
        freedom_case = cases.FreedomCase()
        sup_a = [0, 0, 0]
        sup_c = [0, 0, 0]
        sup_a[0] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=0)
        sup_a[1] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=1)
        sup_a[2] = freedom_case.add_nodal_support(node=self.node_a, val=0, dof=5)
        sup_c[0] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=0)
        sup_c[1] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=1)
        sup_c[2] = freedom_case.add_nodal_support(node=self.node_c, val=0, dof=5)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=self.node_b, val=100e3/np.sqrt(2), dof=0)
        load_case.add_nodal_load(node=self.node_b, val=-100e3/np.sqrt(2), dof=1)
        load_case.add_nodal_load(node=self.node_b, val=50e6, dof=5)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=self.analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_b = self.node_b.get_displacements(analysis_case)
        u_b = disp_b[0]
        v_b = disp_b[1]
        rz_b = disp_b[5]

        self.assertEqual(np.around(u_b, 4), 0.4415)
        self.assertEqual(np.around(v_b, 4), -0.3999)
        self.assertEqual(np.around(rz_b, 5), 0.00169)

    def test_example4_15(self):
        # create materials
        steel = Material(name='steel_imperial', elastic_modulus=29e3, poissons_ratio=0.3, rho=0)

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()
        analysis.post.n_subdiv = 3

        beam = Section(ixx=128.5)

        # create nodes
        n = 20  # number of beam segments (must be an even number to ensure point load in centre)
        L = 300  # length
        k_L = 15  # spring length
        k = 1.5  # foundation modulus
        dx = L/n
        spring = Section(area=k*dx*k_L/29e3)
        beam_nodes = []
        spring_nodes = []

        # beam nodes
        for i in range(n+1):
            beam_nodes.append(analysis.create_node(coords=[i*L/n]))

        # spring nodes
        for i in range(n-1):
            spring_nodes.append(analysis.create_node(coords=[(i+1)*L/n, -k_L]))

        # create elements
        beam_elements = []
        spring_elements = []

        # create beam elements
        for i in range(n):
            beam_elements.append(analysis.create_element(
                el_type='EB2-2D', nodes=[beam_nodes[i], beam_nodes[i+1]], material=steel,
                section=beam
            ))

        # create spring elements
        for i in range(n-1):
            spring_elements.append(analysis.create_element(
                el_type='Bar2-2D', nodes=[beam_nodes[i+1], spring_nodes[i]], material=steel,
                section=spring
            ))

        # add supports
        freedom_case = cases.FreedomCase()
        spring_supports_x = []
        spring_supports_y = []

        # add end supports
        freedom_case.add_nodal_support(node=beam_nodes[0], val=0, dof=0)
        freedom_case.add_nodal_support(node=beam_nodes[0], val=0, dof=1)
        freedom_case.add_nodal_support(node=beam_nodes[-1], val=0, dof=1)

        # add spring support
        for spring_node in spring_nodes:
            spring_supports_x.append(freedom_case.add_nodal_support(
                node=spring_node, val=0, dof=0))
            spring_supports_y.append(freedom_case.add_nodal_support(
                node=spring_node, val=0, dof=1))

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=beam_nodes[int(n/2)], val=-40, dof=1)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp = beam_nodes[int(n/2)].get_displacements(analysis_case)
        v_p = disp[1]

        self.assertEqual(np.around(v_p, 3), -0.238)

        # check bending moments
        (_, m) = beam_elements[int(n/2)].get_bmd(n=2, analysis_case=analysis_case)
        self.assertEqual(np.around(m[0], 0), -547)


if __name__ == "__main__":
    unittest.main()
