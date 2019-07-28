from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.feasolve import SolverSettings

# ------------
# preprocessor
# ------------

# lists
nodes = []  # list holding the node objects
elements = []  # list holding the element objects
freedom_cases = []  # list holding the freedom cases
load_cases = []  # list holding the load cases
analysis_cases = []  # list holding the analysis cases

# create 2d frame analyis object
analysis = FrameAnalysis2D()

# create materials and sections
steel = Steel()
section = Section(area=3230, ixx=23.6e6)

# create nodes
nodes.append(analysis.create_node(coords=[0, 0]))
nodes.append(analysis.create_node(coords=[1500, 0]))
nodes.append(analysis.create_node(coords=[3000, 1500]))

# create beam elements
for i in range(2):
    elements.append(analysis.create_element(
        el_type='EB2-2D',
        nodes=[nodes[i], nodes[i+1]],
        material=steel,
        section=section
    ))

# create 3 freedom cases
for i in range(3):
    freedom_cases.append(cases.FreedomCase())

# freedom case 1: pin middle node, roller right node
freedom_cases[0].add_nodal_support(node=nodes[1], val=0, dof=0)
freedom_cases[0].add_nodal_support(node=nodes[1], val=0, dof=1)
freedom_cases[0].add_nodal_support(node=nodes[2], val=0, dof=1)

# freedom case 2: imposed displacement left node, pin middle node, roller right node
freedom_cases[1].add_nodal_support(node=nodes[0], val=-10, dof=1)
freedom_cases[1].add_nodal_support(node=nodes[1], val=0, dof=0)
freedom_cases[1].add_nodal_support(node=nodes[1], val=0, dof=1)
freedom_cases[1].add_nodal_support(node=nodes[2], val=0, dof=1)

# freedom case 3: imposed displacement and rotation left node, pin middle node, roller right node
freedom_cases[2].add_nodal_support(node=nodes[0], val=0.01, dof=5)
freedom_cases[2].add_nodal_support(node=nodes[1], val=0, dof=0)
freedom_cases[2].add_nodal_support(node=nodes[1], val=0, dof=1)
freedom_cases[2].add_nodal_support(node=nodes[2], val=0, dof=1)

# create 2 load cases
for i in range(2):
    load_cases.append(cases.LoadCase())

# load case 1: point load at left node
load_cases[0].add_nodal_load(node=nodes[0], val=-1e4, dof=1)

# load case 2: no loads

# add analysis cases
analysis_cases.append(cases.AnalysisCase(freedom_case=freedom_cases[0], load_case=load_cases[0]))
analysis_cases.append(cases.AnalysisCase(freedom_case=freedom_cases[1], load_case=load_cases[1]))
analysis_cases.append(cases.AnalysisCase(freedom_case=freedom_cases[2], load_case=load_cases[1]))

# ------
# solver
# ------

settings = SolverSettings()
settings.linear_static.time_info = True

LinearStatic(analysis=analysis, analysis_cases=analysis_cases, solver_settings=settings).solve()

# ----
# post
# ----

# plot the results for each analysis case
for analysis_case in analysis_cases:
    analysis.post.plot_geom(analysis_case=analysis_case)
    analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=10)
    analysis.post.plot_frame_forces(
        analysis_case=analysis_case, axial=True, shear=True, moment=True)
    analysis.post.plot_reactions(analysis_case=analysis_case)
