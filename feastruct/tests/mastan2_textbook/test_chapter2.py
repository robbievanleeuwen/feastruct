import unittest
import numpy as np
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic


class TestChapter2(unittest.TestCase):
    """Tests examples and problems from Chapter 2 of the mastan2 textbook."""

    def setUp(self):
        # create materials
        self.steel = Steel()

    def test_example2_1(self):
        nodes = []  # list holding the node objects
        elements = []  # list holding the element objects

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create sections
        steel = Steel()
        section1 = Section(area=6000)
        section2 = Section(area=8000)

        # create nodes
        nodes.append(analysis.create_node(coords=[0, 0]))
        nodes.append(analysis.create_node(coords=[6000, 4000]))
        nodes.append(analysis.create_node(coords=[9000, 0]))

        # create truss elements
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[0], nodes[1]], material=steel, section=section1
        ))
        elements.append(analysis.create_element(
            el_type='Bar2-2D', nodes=[nodes[1], nodes[2]], material=steel, section=section2
        ))

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=nodes[0], val=0, dof=0)
        freedom_case.add_nodal_support(node=nodes[0], val=0, dof=1)
        freedom_case.add_nodal_support(node=nodes[2], val=0, dof=0)
        freedom_case.add_nodal_support(node=nodes[2], val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=nodes[1], val=500e3, dof=0)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        n1_disp = nodes[1].get_displacements(analysis_case)

        u = n1_disp[0]
        v = n1_disp[1]

        self.assertEqual(np.around(u, 2), 2.41)
        self.assertEqual(np.around(v, 2), 0.72)

        # check axial forces
        (_, n1) = elements[0].get_afd(n=1, analysis_case=analysis_case)
        (_, n2) = elements[1].get_afd(n=1, analysis_case=analysis_case)

        self.assertEqual(np.around(n1/1e3, 1), 400.6)
        self.assertEqual(np.around(n2/1e3, 1), -277.8)

    def test_example2_2(self):
        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create sections
        steel = Steel()
        section1 = Section(area=5000)
        section2 = Section(area=3000)

        # create nodes
        a = analysis.create_node(coords=[0, 0])
        b = analysis.create_node(coords=[3000, 3000])
        c = analysis.create_node(coords=[3000, -3000])
        d = analysis.create_node(coords=[-3000/np.tan(np.pi/6), 3000])
        e = analysis.create_node(coords=[-3000/np.tan(np.pi/6), -3000])

        # create truss elements
        el1 = analysis.create_element(
            el_type='Bar2-2D', nodes=[e, a], material=steel, section=section1)
        el2 = analysis.create_element(
            el_type='Bar2-2D', nodes=[d, a], material=steel, section=section1)
        el3 = analysis.create_element(
            el_type='Bar2-2D', nodes=[b, a], material=steel, section=section2)
        el4 = analysis.create_element(
            el_type='Bar2-2D', nodes=[c, a], material=steel, section=section2)

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=b, val=0, dof=0)
        freedom_case.add_nodal_support(node=b, val=0, dof=1)
        freedom_case.add_nodal_support(node=c, val=0, dof=0)
        freedom_case.add_nodal_support(node=c, val=0, dof=1)
        freedom_case.add_nodal_support(node=d, val=0, dof=0)
        freedom_case.add_nodal_support(node=d, val=0, dof=1)
        freedom_case.add_nodal_support(node=e, val=0, dof=0)
        freedom_case.add_nodal_support(node=e, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=a, val=1000e3, dof=0)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check node displacement
        disp_a = a.get_displacements(analysis_case)
        u = disp_a[0]
        v = disp_a[1]

        self.assertEqual(np.around(u, 2), 2.55)
        self.assertEqual(np.around(v, 2), 0)

        # check axial forces
        (_, n1) = el1.get_afd(n=1, analysis_case=analysis_case)
        (_, n2) = el2.get_afd(n=1, analysis_case=analysis_case)
        (_, n3) = el3.get_afd(n=1, analysis_case=analysis_case)
        (_, n4) = el4.get_afd(n=1, analysis_case=analysis_case)

        self.assertEqual(np.around(n1/1e3, 1), 368.8)
        self.assertEqual(np.around(n2/1e3, 1), 368.8)
        self.assertEqual(np.around(n3/1e3, 1), -255.5)
        self.assertEqual(np.around(n4/1e3, 1), -255.5)

    def test_example2_3(self):
        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create sections
        steel = Steel()
        section1 = Section(area=4000)
        section2 = Section(area=1000)

        # create nodes
        a = analysis.create_node(coords=[0, 0])
        b = analysis.create_node(coords=[2000, 2000])
        c = analysis.create_node(coords=[-4000, 2000])

        # create truss elements
        el1 = analysis.create_element(
            el_type='Bar2-2D', nodes=[c, a], material=steel, section=section1)
        el2 = analysis.create_element(
            el_type='Bar2-2D', nodes=[a, b], material=steel, section=section2)

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=b, val=0, dof=0)
        freedom_case.add_nodal_support(node=b, val=0, dof=1)
        freedom_case.add_nodal_support(node=c, val=0, dof=0)
        freedom_case.add_nodal_support(node=c, val=0, dof=1)
        sup1 = freedom_case.add_nodal_support(node=a, val=0, dof=0)
        sup2 = freedom_case.add_nodal_support(node=a, val=-5, dof=1)

        # add loads
        load_case = cases.LoadCase()

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check reactions
        self.assertEqual(np.around(sup1.get_reaction(analysis_case)/1e3, 1), 181.0)
        self.assertEqual(np.around(sup2.get_reaction(analysis_case)/1e3, 1), -355.7)


if __name__ == "__main__":
    unittest.main()
