import unittest
import numpy as np
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic


class TestChapter3(unittest.TestCase):
    """Tests examples and problems from Chapter 3 of the mastan2 textbook."""

    def setUp(self):
        # create materials
        self.steel = Steel()

    def test_example3_2(self):
        nodes = []  # list holding the node objects
        elements = []  # list holding the element objects

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create sections
        steel = Steel()
        section1 = Section(area=15e3)
        section2 = Section(area=18e3)
        section3 = Section(area=20e3)

        # create nodes
        nodes.append(analysis.create_node(coords=[0, 0]))
        nodes.append(analysis.create_node(coords=[4000/np.tan(np.pi/6), 4000]))
        nodes.append(analysis.create_node(coords=[4000/np.tan(np.pi/6) + 4000, 0]))

        # create truss elements
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[0], nodes[1]], material=steel, section=section1
        ))
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[1], nodes[2]], material=steel, section=section3
        ))
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[0], nodes[2]], material=steel, section=section2
        ))

        # add supports
        freedom_case = cases.FreedomCase()
        sup1 = freedom_case.add_nodal_support(node=nodes[0], val=0, dof=0)
        sup2 = freedom_case.add_nodal_support(node=nodes[0], val=0, dof=1)
        sup3 = freedom_case.add_nodal_support(node=nodes[2], val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=nodes[1], val=500e3*np.cos(2*np.pi/9), dof=0)
        load_case.add_nodal_load(node=nodes[1], val=500e3*np.sin(2*np.pi/9), dof=1)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_a = nodes[1].get_displacements(analysis_case)
        disp_b = nodes[2].get_displacements(analysis_case)

        u_a = disp_a[0]
        v_a = disp_a[1]
        u_b = disp_b[0]

        self.assertEqual(np.around(u_a, 2), 0.87)
        self.assertEqual(np.around(v_a, 2), 1.24)
        self.assertEqual(np.around(u_b, 2), -0.19)

        # check axial forces
        (_, n1) = elements[0].get_afd(n=1, analysis_case=analysis_case)
        (_, n2) = elements[1].get_afd(n=1, analysis_case=analysis_case)
        (_, n3) = elements[2].get_afd(n=1, analysis_case=analysis_case)

        self.assertEqual(np.around(n1/1e3, 1), 515.7)
        self.assertEqual(np.around(n2/1e3, 1), 89.9)
        self.assertEqual(np.around(n3/1e3, 1), -63.6)

        # check reactions
        rc_x = sup1.get_reaction(analysis_case=analysis_case)
        rc_y = sup2.get_reaction(analysis_case=analysis_case)
        rb_y = sup3.get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(rc_x/1e3, 1), -383.0)
        self.assertEqual(np.around(rc_y/1e3, 1), -257.8)
        self.assertEqual(np.around(rb_y/1e3, 1), -63.6)

    def test_example3_3(self):
        nodes = []  # list holding the node objects
        elements = []  # list holding the element objects

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create sections
        steel = Steel()
        section1 = Section(area=15e3)
        section2 = Section(area=18e3)
        section3 = Section(area=20e3)

        # create nodes
        nodes.append(analysis.create_node(coords=[0, 0]))
        nodes.append(analysis.create_node(coords=[4000/np.tan(np.pi/6), 4000]))
        nodes.append(analysis.create_node(coords=[4000/np.tan(np.pi/6) + 4000, 0]))
        nodes.append(analysis.create_node(coords=[0, 4000]))

        # create truss elements
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[0], nodes[1]], material=steel, section=section1
        ))
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[1], nodes[2]], material=steel, section=section3
        ))
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[0], nodes[2]], material=steel, section=section2
        ))
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[3], nodes[1]], material=steel, section=section3
        ))

        # add supports
        freedom_case = cases.FreedomCase()
        sup1 = freedom_case.add_nodal_support(node=nodes[0], val=0, dof=0)
        sup2 = freedom_case.add_nodal_support(node=nodes[0], val=0, dof=1)
        sup3 = freedom_case.add_nodal_support(node=nodes[2], val=0, dof=1)
        sup4 = freedom_case.add_nodal_support(node=nodes[3], val=0, dof=0)
        sup5 = freedom_case.add_nodal_support(node=nodes[3], val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=nodes[1], val=500e3*np.cos(2*np.pi/9), dof=0)
        load_case.add_nodal_load(node=nodes[1], val=500e3*np.sin(2*np.pi/9), dof=1)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_a = nodes[1].get_displacements(analysis_case)
        disp_b = nodes[2].get_displacements(analysis_case)

        u_a = disp_a[0]
        v_a = disp_a[1]
        u_b = disp_b[0]

        self.assertEqual(np.around(u_a, 2), 0.38)
        self.assertEqual(np.around(v_a, 2), 1.23)
        self.assertEqual(np.around(u_b, 2), -0.44)

        # check reactions
        rc_x = sup1.get_reaction(analysis_case=analysis_case)
        rc_y = sup2.get_reaction(analysis_case=analysis_case)
        rb_y = sup3.get_reaction(analysis_case=analysis_case)
        rd_x = sup4.get_reaction(analysis_case=analysis_case)
        rd_y = sup5.get_reaction(analysis_case=analysis_case)

        self.assertEqual(np.around(rc_x/1e3, 1), -162.5)
        self.assertEqual(np.around(rc_y/1e3, 1), -177.1)
        self.assertEqual(np.around(rb_y/1e3, 1), -144.3)
        self.assertEqual(np.around(rd_x/1e3, 1), -220.5)
        self.assertEqual(np.around(rd_y/1e3, 1), 0)


if __name__ == "__main__":
    unittest.main()
