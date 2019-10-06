from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.feasolve import SolverSettings

# ------------
# preprocessor
# ------------

# constants & lists
length = 5000  # length of the beam
udl = -10  # value of the udl

# create 2d frame analysis object
analysis = FrameAnalysis2D()

# create materials and sections
steel = Steel()
section = Section(area=3230, ixx=23.6e6)

# create nodes
node1 = analysis.create_node(coords=[0])
node2 = analysis.create_node(coords=[length])
node3 = analysis.create_node(coords=[2*length])

# create beam element
beam1 = analysis.create_element(
    el_type='EB2-2D', nodes=[node1, node2], material=steel, section=section
)
beam2 = analysis.create_element(
    el_type='EB2-2D', nodes=[node2, node3], material=steel, section=section
)

# add supports
freedom_case = cases.FreedomCase()
freedom_case.add_nodal_support(node=node1, val=0, dof=0)
freedom_case.add_nodal_support(node=node1, val=0, dof=1)
freedom_case.add_nodal_support(node=node2, val=0, dof=1)
freedom_case.add_nodal_support(node=node3, val=0, dof=1)

# add loads
load_case = cases.LoadCase()
load_case.add_element_load(beam1.generate_udl(q=udl))
load_case.add_element_load(beam2.generate_udl(q=udl))

# add analysis case
analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

# ------
# solver
# ------

settings = SolverSettings()
settings.linear_static.time_info = True

LinearStatic(analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings).solve()

# ----
# post
# ----

analysis.post.plot_geom(analysis_case=analysis_case)
analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=100)
analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)
analysis.post.plot_reactions(analysis_case=analysis_case)
