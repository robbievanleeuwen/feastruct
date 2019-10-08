import unittest
import numpy as np
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic


class TestUDL(unittest.TestCase):
    """Tests problems related to 1D beam bending from the American Wood Council:
    https://www.awc.org/pdf/codes-standards/publications/design-aids/AWC-DA6-BeamFormulas-0710.pdf
    """

    def setUp(self):
        self.steel = Steel()
        self.elastic_modulus = self.steel.elastic_modulus
        self.ixx = np.random.uniform(10e6, 200e6)
        self.length = np.random.uniform(2e3, 10e3)
        self.q = -np.random.uniform(1, 10)
        self.pl = -np.random.uniform(5e3, 50e3)

    def test_fig1(self):
        """Simple Beam – Uniformly Distributed Load"""

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create section
        section = Section(ixx=self.ixx)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[self.length])

        # create beam elements
        element = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=self.steel, section=section
        )

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=node_a, val=0, dof=0)
        freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        freedom_case.add_nodal_support(node=node_b, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_element_load(element.generate_udl(q=self.q))

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check displacements
        def analytical_disp(x):
            factor = self.q * x / 24 / self.elastic_modulus / self.ixx
            l0 = self.length

            return factor * (l0 * l0 * l0 - 2 * l0 * x * x + x * x * x)

        # get displacements
        displacements = element.get_displacements(11, analysis_case)

        # loop through each station
        for disp in displacements:
            xi = disp[0]
            x = self.length * xi
            v = disp[2]

            # check displacements
            self.assertTrue(np.isclose(v, analytical_disp(x), atol=1e-06))

        # check max displacement
        l0 = self.length
        v_max = 5 * self.q * l0 * l0 * l0 * l0 / 384 / self.elastic_modulus / self.ixx

        # check value
        self.assertTrue(np.isclose(abs(v_max), max(np.abs(displacements[:, 2]))))

        # check position
        self.assertTrue(np.isclose(0.5, displacements[np.abs(displacements[:, 2]).argmax(), 0],
                                   atol=1e-06))

        # check bending moments
        def analytical_bmd(x):
            return self.q * x / 2 * (self.length - x)

        # get bmd
        (xis, bmd) = element.get_bmd(11, analysis_case)

        # loop through each station
        for (i, m) in enumerate(bmd):
            xi = xis[i]
            x = self.length * xi

            # check bending moment
            self.assertTrue(np.isclose(m, analytical_bmd(x), atol=1e-06))

        # check max bending moment
        l0 = self.length
        m_max = self.q * l0 * l0 / 8

        # check value
        self.assertTrue(np.isclose(abs(m_max), max(np.abs(bmd)), atol=1e-06))

        # check position
        self.assertTrue(np.isclose(0.5, xis[np.abs(bmd).argmax()], atol=1e-06))

        # check shear force
        def analytical_sfd(x):
            return self.q * (x - self.length / 2)

        # get sfd
        (xis, sfd) = element.get_sfd(11, analysis_case)

        # loop through each station
        for (i, sf) in enumerate(sfd):
            xi = xis[i]
            x = self.length * xi

            # check shear force
            self.assertTrue(np.isclose(sf, analytical_sfd(x), atol=1e-06))

    def test_fig2(self):
        """Simple Beam – Uniform Load Partially Distributed"""

        a = self.length * np.random.uniform(0.1, 0.4)
        c = self.length * np.random.uniform(0.1, 0.4)
        b = self.length - a - c

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create section
        section = Section(ixx=self.ixx)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[a])
        node_c = analysis.create_node(coords=[a+b])
        node_d = analysis.create_node(coords=[self.length])

        # create beam elements
        element_ab = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=self.steel, section=section
        )
        element_bc = analysis.create_element(
            el_type='EB2-2D', nodes=[node_b, node_c], material=self.steel, section=section
        )
        element_cd = analysis.create_element(
            el_type='EB2-2D', nodes=[node_c, node_d], material=self.steel, section=section
        )

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=node_a, val=0, dof=0)
        sup1 = freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        sup2 = freedom_case.add_nodal_support(node=node_d, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_element_load(element_bc.generate_udl(q=self.q))

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check reactions
        r1 = -sup1.get_reaction(analysis_case)
        r2 = -sup2.get_reaction(analysis_case)

        self.assertTrue(np.isclose(r1, self.q * b / 2 / self.length * (2 * c + b), atol=1e-06))
        self.assertTrue(np.isclose(r2, self.q * b / 2 / self.length * (2 * a + b), atol=1e-06))

        # check bending moments
        def analytical_bmd_ab(x):
            return r1 * x

        def analytical_bmd_bc(x):
            return r1 * x - self.q / 2 * (x - a) * (x - a)

        def analytical_bmd_cd(x):
            return r2 * (self.length - x)

        # get bmds
        (xis_ab, bmd_ab) = element_ab.get_bmd(11, analysis_case)
        (xis_bc, bmd_bc) = element_bc.get_bmd(11, analysis_case)
        (xis_cd, bmd_cd) = element_cd.get_bmd(11, analysis_case)

        # element_ab - loop through each station
        for (i, m) in enumerate(bmd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_ab(x), atol=1e-06))

        # element_bc - loop through each station
        for (i, m) in enumerate(bmd_bc):
            xi = xis_bc[i]
            x = b * xi + a

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_bc(x), atol=1e-06))

        # element_cd - loop through each station
        for (i, m) in enumerate(bmd_cd):
            xi = xis_cd[i]
            x = c * xi + a + b

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_cd(x), atol=1e-06))

        # check max bending moment
        m_max = r1 * (a + r1 / 2 / self.q)
        pos = a + r1 / self.q
        x = 1 / b * (pos - a)

        # check value
        self.assertTrue(np.isclose(abs(m_max), max(np.abs(bmd_bc)), atol=1e-06))

        # check position
        self.assertTrue(np.isclose(x, xis_bc[np.abs(bmd_bc).argmax()], atol=1e-06))

        # check shear force
        def analytical_sfd_ab(x):
            return -r1

        def analytical_sfd_bc(x):
            return -r1 + self.q * (x - a)

        def analytical_sfd_cd(x):
            return r2

        # get sfds
        (xis_ab, sfd_ab) = element_ab.get_sfd(11, analysis_case)
        (xis_bc, sfd_bc) = element_bc.get_sfd(11, analysis_case)
        (xis_cd, sfd_cd) = element_cd.get_sfd(11, analysis_case)

        # element_ab - loop through each station
        for (i, sf) in enumerate(sfd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_ab(x), atol=1e-06))

        # element_bc - loop through each station
        for (i, sf) in enumerate(sfd_bc):
            xi = xis_bc[i]
            x = b * xi + a

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_bc(x), atol=1e-06))

        # element_cd - loop through each station
        for (i, sf) in enumerate(sfd_cd):
            xi = xis_cd[i]
            x = c * xi + a + b

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_cd(x), atol=1e-06))

    def test_fig3(self):
        """Simple Beam – Uniform Load Partially Distributed at One End"""

        a = self.length * np.random.uniform(0.1, 0.9)

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create section
        section = Section(ixx=self.ixx)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[a])
        node_c = analysis.create_node(coords=[self.length])

        # create beam elements
        element_ab = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=self.steel, section=section
        )
        element_bc = analysis.create_element(
            el_type='EB2-2D', nodes=[node_b, node_c], material=self.steel, section=section
        )

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=node_a, val=0, dof=0)
        sup1 = freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        sup2 = freedom_case.add_nodal_support(node=node_c, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_element_load(element_ab.generate_udl(q=self.q))

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check reactions
        r1 = -sup1.get_reaction(analysis_case)
        r2 = -sup2.get_reaction(analysis_case)

        self.assertTrue(np.isclose(r1, self.q * a / 2 / self.length * (2 * self.length - a),
                                   atol=1e-06))
        self.assertTrue(np.isclose(r2, self.q * a * a / 2 / self.length, atol=1e-06))

        # check displacements
        def analytical_disp_ab(x):
            l0 = self.length
            factor = self.q * x / 24 / self.elastic_modulus / self.ixx / l0

            return factor * (a * a * (2 * l0 - a) * (2 * l0 - a) - 2 * a * x * x * (
                2 * l0 - a) + l0 * x * x * x)

        def analytical_disp_bc(x):
            l0 = self.length
            factor = self.q * a * a * (l0 - x) / 24 / self.elastic_modulus / self.ixx / l0

            return factor * (4 * x * l0 - 2 * x * x - a * a)

        # get displacements
        displacements_ab = element_ab.get_displacements(11, analysis_case)
        displacements_bc = element_bc.get_displacements(11, analysis_case)

        # loop through each station
        for disp in displacements_ab:
            xi = disp[0]
            x = a * xi
            v = disp[2]

            # check displacements
            self.assertTrue(np.isclose(v, analytical_disp_ab(x), atol=1e-06))

        # loop through each station
        for disp in displacements_bc:
            xi = disp[0]
            x = (self.length - a) * xi + a
            v = disp[2]

            # check displacements
            self.assertTrue(np.isclose(v, analytical_disp_bc(x), atol=1e-06))

        # check bending moments
        def analytical_bmd_ab(x):
            return r1 * x - self.q * x * x / 2

        def analytical_bmd_bc(x):
            return r2 * (self.length - x)

        # get bmds
        (xis_ab, bmd_ab) = element_ab.get_bmd(11, analysis_case)
        (xis_bc, bmd_bc) = element_bc.get_bmd(11, analysis_case)

        # element_ab - loop through each station
        for (i, m) in enumerate(bmd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_ab(x), atol=1e-06))

        # element_bc - loop through each station
        for (i, m) in enumerate(bmd_bc):
            xi = xis_bc[i]
            x = (self.length - a) * xi + a

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_bc(x), atol=1e-06))

        # check max bending moment
        m_max = r1 * r1 / 2 / self.q
        pos = r1 / self.q
        x = pos / a

        # check value
        self.assertTrue(np.isclose(abs(m_max), max(np.abs(bmd_ab)), atol=1e-06))

        # check position
        self.assertTrue(np.isclose(x, xis_ab[np.abs(bmd_ab).argmax()], atol=1e-06))

        # check shear force
        def analytical_sfd_ab(x):
            return -r1 + self.q * x

        def analytical_sfd_bc(x):
            return r2

        # get sfds
        (xis_ab, sfd_ab) = element_ab.get_sfd(11, analysis_case)
        (xis_bc, sfd_bc) = element_bc.get_sfd(11, analysis_case)

        # element_ab - loop through each station
        for (i, sf) in enumerate(sfd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_ab(x), atol=1e-06))

        # element_bc - loop through each station
        for (i, sf) in enumerate(sfd_bc):
            xi = xis_bc[i]
            x = (self.length - a) * xi + a

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_bc(x), atol=1e-06))

    def test_fig4(self):
        """Simple Beam – Uniform Load Partially Distributed at Each End"""

        a = self.length * np.random.uniform(0.1, 0.4)
        c = self.length * np.random.uniform(0.1, 0.4)
        b = self.length - a - c
        q2 = -np.random.uniform(1, 10)

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create section
        section = Section(ixx=self.ixx)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[a])
        node_c = analysis.create_node(coords=[a+b])
        node_d = analysis.create_node(coords=[self.length])

        # create beam elements
        element_ab = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=self.steel, section=section
        )
        element_bc = analysis.create_element(
            el_type='EB2-2D', nodes=[node_b, node_c], material=self.steel, section=section
        )
        element_cd = analysis.create_element(
            el_type='EB2-2D', nodes=[node_c, node_d], material=self.steel, section=section
        )

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=node_a, val=0, dof=0)
        sup1 = freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        sup2 = freedom_case.add_nodal_support(node=node_d, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_element_load(element_ab.generate_udl(q=self.q))
        load_case.add_element_load(element_cd.generate_udl(q=q2))

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check reactions
        r1 = -sup1.get_reaction(analysis_case)
        r1_ana = (self.q * a * (2 * self.length - a) + q2 * c * c) / (2 * self.length)
        r2 = -sup2.get_reaction(analysis_case)
        r2_ana = (q2 * c * (2 * self.length - c) + self.q * a * a) / (2 * self.length)

        self.assertTrue(np.isclose(r1, r1_ana, atol=1e-06))
        self.assertTrue(np.isclose(r2, r2_ana, atol=1e-06))

        # check bending moments
        def analytical_bmd_ab(x):
            return r1 * x - self.q * 0.5 * x * x

        def analytical_bmd_bc(x):
            return r1 * x - self.q * a * 0.5 * (2 * x - a)

        def analytical_bmd_cd(x):
            return r2 * (self.length - x) - q2 * (self.length - x) * (self.length - x) * 0.5

        # get bmds
        (xis_ab, bmd_ab) = element_ab.get_bmd(11, analysis_case)
        (xis_bc, bmd_bc) = element_bc.get_bmd(11, analysis_case)
        (xis_cd, bmd_cd) = element_cd.get_bmd(11, analysis_case)

        # element_ab - loop through each station
        for (i, m) in enumerate(bmd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_ab(x), atol=1e-06))

        # element_bc - loop through each station
        for (i, m) in enumerate(bmd_bc):
            xi = xis_bc[i]
            x = b * xi + a

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_bc(x), atol=1e-06))

        # element_cd - loop through each station
        for (i, m) in enumerate(bmd_cd):
            xi = xis_cd[i]
            x = c * xi + a + b

            # check bending moments
            self.assertTrue(np.isclose(m, analytical_bmd_cd(x), atol=1e-06))

        # check max bending moment
        if abs(r1) < abs(self.q * a):
            m_max = r1 * r1 / 2 / self.q
            pos = r1 / self.q
            x = pos / a

            # check value
            self.assertTrue(np.isclose(abs(m_max), max(np.abs(bmd_ab)), atol=1e-06))

            # check position
            self.assertTrue(np.isclose(x, xis_ab[np.abs(bmd_ab).argmax()], atol=1e-06))

        if abs(r2) < abs(q2 * c):
            m_max = r2 * r2 / 2 / q2
            pos = self.length - r2 / q2
            x = 1 / c * (pos - a - b)

            # check value
            self.assertTrue(np.isclose(abs(m_max), max(np.abs(bmd_cd)), atol=1e-06))

            # check position
            self.assertTrue(np.isclose(x, xis_cd[np.abs(bmd_cd).argmax()], atol=1e-06))

        # check shear force
        def analytical_sfd_ab(x):
            return -r1 + self.q * x

        def analytical_sfd_bc(x):
            return -r1 + self.q * a

        def analytical_sfd_cd(x):
            return r2 - q2 * (self.length - x)

        # get sfds
        (xis_ab, sfd_ab) = element_ab.get_sfd(11, analysis_case)
        (xis_bc, sfd_bc) = element_bc.get_sfd(11, analysis_case)
        (xis_cd, sfd_cd) = element_cd.get_sfd(11, analysis_case)

        # element_ab - loop through each station
        for (i, sf) in enumerate(sfd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_ab(x), atol=1e-06))

        # element_bc - loop through each station
        for (i, sf) in enumerate(sfd_bc):
            xi = xis_bc[i]
            x = b * xi + a

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_bc(x), atol=1e-06))

        # element_cd - loop through each station
        for (i, sf) in enumerate(sfd_cd):
            xi = xis_cd[i]
            x = c * xi + a + b

            # check shear forces
            self.assertTrue(np.isclose(sf, analytical_sfd_cd(x), atol=1e-06))

    def test_fig5(self):
        """Simple Beam – Load Increasing Uniformly to One End"""

        # not yet implemented
        pass

    def test_fig6(self):
        """Simple Beam – Load Increasing Uniformly to Center"""

        # not yet implemented
        pass

    def test_fig7(self):
        """Simple Beam – Concentrated Load at Center"""

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create section
        section = Section(ixx=self.ixx)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[self.length * 0.5])
        node_c = analysis.create_node(coords=[self.length])

        # create beam elements
        element_ab = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=self.steel, section=section
        )
        element_bc = analysis.create_element(
            el_type='EB2-2D', nodes=[node_b, node_c], material=self.steel, section=section
        )

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=node_a, val=0, dof=0)
        freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        freedom_case.add_nodal_support(node=node_c, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=node_b, val=self.pl, dof=1)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check displacements
        def analytical_disp_ab(x):
            factor = self.pl * x / 48 / self.elastic_modulus / self.ixx
            l0 = self.length

            return factor * (3 * l0 * l0 - 4 * x * x)

        def analytical_disp_bc(x):
            x = self.length - x

            factor = self.pl * x / 48 / self.elastic_modulus / self.ixx
            l0 = self.length

            return factor * (3 * l0 * l0 - 4 * x * x)

        # get displacements
        displacements_ab = element_ab.get_displacements(11, analysis_case)
        displacements_bc = element_bc.get_displacements(11, analysis_case)

        # loop through each station
        for disp in displacements_ab:
            xi = disp[0]
            x = self.length * 0.5 * xi
            v = disp[2]

            # check displacements
            self.assertTrue(np.isclose(v, analytical_disp_ab(x), atol=1e-06))

        # loop through each station
        for disp in displacements_bc:
            xi = disp[0]
            x = self.length * 0.5 + self.length * 0.5 * xi
            v = disp[2]

            # check displacements
            self.assertTrue(np.isclose(v, analytical_disp_bc(x), atol=1e-06))

        # check max displacement
        l0 = self.length
        v_max = self.pl * l0 * l0 * l0 / 48 / self.elastic_modulus / self.ixx

        # check value
        self.assertTrue(np.isclose(abs(v_max), max(np.abs(displacements_ab[:, 2])), atol=1e-06))

        # check position
        self.assertTrue(
            np.isclose(1, displacements_ab[np.abs(displacements_ab[:, 2]).argmax(), 0],
                       atol=1e-06))

        # check bending moments
        def analytical_bmd_ab(x):
            return self.pl * x / 2

        def analytical_bmd_bc(x):
            x = self.length - x
            return self.pl * x / 2

        # get bmd
        (xis_ab, bmd_ab) = element_ab.get_bmd(11, analysis_case)
        (xis_bc, bmd_bc) = element_bc.get_bmd(11, analysis_case)

        # loop through each station
        for (i, m) in enumerate(bmd_ab):
            xi = xis_ab[i]
            x = self.length * 0.5 * xi

            # check bending moment
            self.assertTrue(np.isclose(m, analytical_bmd_ab(x), atol=1e-06))

        # loop through each station
        for (i, m) in enumerate(bmd_bc):
            xi = xis_bc[i]
            x = self.length * 0.5 + self.length * 0.5 * xi

            # check bending moment
            self.assertTrue(np.isclose(m, analytical_bmd_bc(x), atol=1e-06))

        # check max bending moment
        l0 = self.length
        m_max = self.pl * l0 / 4

        # check value
        self.assertTrue(np.isclose(abs(m_max), max(np.abs(bmd_ab)), atol=1e-06))

        # check position
        self.assertTrue(np.isclose(1, xis_ab[np.abs(bmd_ab).argmax()], atol=1e-06))

        # check shear force
        def analytical_sfd_ab(x):
            return -self.pl / 2

        def analytical_sfd_bc(x):
            return self.pl / 2

        # get sfd
        (xis, sfd_ab) = element_ab.get_sfd(11, analysis_case)
        (xis, sfd_bc) = element_bc.get_sfd(11, analysis_case)

        # loop through each station
        for (i, sf) in enumerate(sfd_ab):
            xi = xis_ab[i]
            x = self.length * 0.5 * xi

            # check shear force
            self.assertTrue(np.isclose(sf, analytical_sfd_ab(x), atol=1e-06))

        # loop through each station
        for (i, sf) in enumerate(sfd_bc):
            xi = xis_bc[i]
            x = self.length * 0.5 + self.length * 0.5 * xi

            # check shear force
            self.assertTrue(np.isclose(sf, analytical_sfd_bc(x), atol=1e-06))

    def test_fig8(self):
        """Simple Beam – Concentrated Load at Any Point"""

        a = self.length * np.random.uniform(0.1, 0.9)
        b = self.length - a

        # create 2d frame analysis object
        analysis = FrameAnalysis2D()

        # create section
        section = Section(ixx=self.ixx)

        # create nodes
        node_a = analysis.create_node(coords=[0])
        node_b = analysis.create_node(coords=[a])
        node_c = analysis.create_node(coords=[self.length])

        # create beam elements
        element_ab = analysis.create_element(
            el_type='EB2-2D', nodes=[node_a, node_b], material=self.steel, section=section
        )
        element_bc = analysis.create_element(
            el_type='EB2-2D', nodes=[node_b, node_c], material=self.steel, section=section
        )

        # add supports
        freedom_case = cases.FreedomCase()
        freedom_case.add_nodal_support(node=node_a, val=0, dof=0)
        freedom_case.add_nodal_support(node=node_a, val=0, dof=1)
        freedom_case.add_nodal_support(node=node_c, val=0, dof=1)

        # add loads
        load_case = cases.LoadCase()
        load_case.add_nodal_load(node=node_b, val=self.pl, dof=1)

        # add analysis case
        analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

        # linear static solver
        LinearStatic(analysis=analysis, analysis_cases=[analysis_case]).solve()

        # check displacements
        def analytical_disp_ab(x):
            factor = self.pl * b * x / 6 / self.elastic_modulus / self.ixx / self.length
            l0 = self.length

            return factor * (l0 * l0 - b * b - x * x)

        def analytical_disp_bc(x):
            l0 = self.length
            factor = self.pl * a * (l0 - x) / 6 / self.elastic_modulus / self.ixx / l0

            return factor * (2 * l0 * x - a * a - x * x)

        # get displacements
        displacements_ab = element_ab.get_displacements(21, analysis_case)
        displacements_bc = element_bc.get_displacements(21, analysis_case)

        # loop through each station
        for disp in displacements_ab:
            xi = disp[0]
            x = a * xi
            v = disp[2]

            # check displacements
            self.assertTrue(np.isclose(v, analytical_disp_ab(x), atol=1e-06))

        # loop through each station
        for disp in displacements_bc:
            xi = disp[0]
            x = b * xi + a
            v = disp[2]

            # check displacements
            self.assertTrue(np.isclose(v, analytical_disp_bc(x), atol=1e-06))

        # check max displacement
        if a > b:
            aa = a
            bb = b
            disps = displacements_ab
        else:
            aa = b
            bb = a
            disps = displacements_bc

        l0 = self.length
        v_max = (self.pl * aa * bb * (aa + 2 * bb) * np.sqrt(3 * aa * (aa + 2 * bb))) / (
            27 * self.elastic_modulus * self.ixx * l0)

        # check value
        self.assertTrue(np.isclose(abs(v_max), max(np.abs(disps[:, 2])), atol=1e-06))

        # check position of max displacement
        pos = np.sqrt((aa * (aa + 2 * bb)) / 3)
        if b > a:
            pos = self.length - pos
            x = 1 / b * (pos - a)
        else:
            x = pos / a

        self.assertTrue(np.isclose(x, disps[np.abs(disps[:, 2]).argmax(), 0], atol=1e-06))

        # check displacement at point load
        v_ab = self.pl * a * a * b * b / 3 / self.elastic_modulus / self.ixx / self.length
        self.assertTrue(np.isclose(v_ab, displacements_bc[0, 2], atol=1e-06))

        # check bending moments
        def analytical_bmd_ab(x):
            return self.pl * b * x / self.length

        def analytical_bmd_bc(x):
            x = self.length - x
            return self.pl * a * x / self.length

        # get bmd
        (xis_ab, bmd_ab) = element_ab.get_bmd(11, analysis_case)
        (xis_bc, bmd_bc) = element_bc.get_bmd(11, analysis_case)

        # loop through each station
        for (i, m) in enumerate(bmd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check bending moment
            self.assertTrue(np.isclose(m, analytical_bmd_ab(x), atol=1e-06))

        # loop through each station
        for (i, m) in enumerate(bmd_bc):
            xi = xis_bc[i]
            x = a + b * xi

            # check bending moment
            self.assertTrue(np.isclose(m, analytical_bmd_bc(x), atol=1e-06))

        # check max bending moment
        m_max = self.pl * a * b / self.length

        # check value
        self.assertTrue(np.isclose(abs(m_max), max(np.abs(bmd_ab)), atol=1e-06))

        # check position
        self.assertTrue(np.isclose(1, xis_ab[np.abs(bmd_ab).argmax()], atol=1e-06))

        # check shear force
        def analytical_sfd_ab(x):
            return -self.pl * b / self.length

        def analytical_sfd_bc(x):
            return self.pl * a / self.length

        # get sfd
        (xis, sfd_ab) = element_ab.get_sfd(11, analysis_case)
        (xis, sfd_bc) = element_bc.get_sfd(11, analysis_case)

        # loop through each station
        for (i, sf) in enumerate(sfd_ab):
            xi = xis_ab[i]
            x = a * xi

            # check shear force
            self.assertTrue(np.isclose(sf, analytical_sfd_ab(x), atol=1e-06))

        # loop through each station
        for (i, sf) in enumerate(sfd_bc):
            xi = xis_bc[i]
            x = a + b * xi

            # check shear force
            self.assertTrue(np.isclose(sf, analytical_sfd_bc(x), atol=1e-06))

    def test_fig9(self):
        """Simple Beam – Two Equal Concentrated Loads Symmetrically Placed"""

        pass


if __name__ == "__main__":
    unittest.main()
