from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame import FrameAnalysis2D
from feastruct.solvers.naturalfrequency import NaturalFrequency
from feastruct.solvers.feasolve import SolverSettings

# ------------
# preprocessor
# ------------

# constants & lists
height = 20000  # length of the beam
num_nodes = 21  # number of nodes to use
num_modes = 6  # number of modes to calculate
nodes = []  # list holding the node objects
elements = []  # list holding the element objects

# create 2d frame analyis object
analysis = FrameAnalysis2D()

# create materials and sections
steel = Steel()
section = Section(area=10400, ixx=88.3e6)

# create nodes
for i in range(num_nodes):
    nodes.append(analysis.create_node(coords=[0, height / (num_nodes - 1) * i]))

# create beam elements
for i in range(num_nodes - 1):
    elements.append(analysis.create_element(
        el_type='EB2-2D',
        nodes=[nodes[i], nodes[i+1]],
        material=steel,
        section=section
    ))

# add supports - fixed base
freedom_case = cases.FreedomCase()
freedom_case.add_nodal_support(node=nodes[0], val=0, dof=0)
freedom_case.add_nodal_support(node=nodes[0], val=0, dof=1)
freedom_case.add_nodal_support(node=nodes[0], val=0, dof=5)

# add analysis case
analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=cases.LoadCase())

# ----------------
# frequency solver
# ----------------

settings = SolverSettings()
settings.natural_frequency.time_info = True
settings.natural_frequency.num_modes = num_modes

solver = NaturalFrequency(
    analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings)
solver.solve()

# --------------
# frequency post
# --------------

for i in range(num_modes):
    analysis.post.plot_frequency_results(analysis_case=analysis_case, frequency_mode=i)
