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
num_nodes = 3  # number of nodes to use
nodes = []  # list holding the node objects
elements = []  # list holding the element objects

# create 2d frame analysis object
analysis = FrameAnalysis2D()

# create materials and sections
steel = Steel()
section = Section(area=3230, ixx=23.6e6)

# create nodes
for i in range(num_nodes):
    nodes.append(analysis.create_node(coords=[length / (num_nodes - 1) * i]))

# create beam elements
for i in range(num_nodes - 1):
    elements.append(analysis.create_element(
        el_type='EB2-2D',
        nodes=[nodes[i], nodes[i+1]],
        material=steel,
        section=section
    ))

# add supports
freedom_case = cases.FreedomCase()
freedom_case.add_nodal_support(node=nodes[0], val=0, dof=0)
freedom_case.add_nodal_support(node=nodes[0], val=0, dof=1)
freedom_case.add_nodal_support(node=nodes[-1], val=0, dof=1)

# add loads
load_case = cases.LoadCase()

for i in range(num_nodes - 2):
    load_case.add_nodal_load(node=nodes[i+1], val=-1e4, dof=1)

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
analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=25)
analysis.post.plot_frame_forces(analysis_case=analysis_case, axial=True)
analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)
analysis.post.plot_reactions(analysis_case=analysis_case)
