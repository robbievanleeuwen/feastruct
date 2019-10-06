import numpy as np
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.linbuckling import LinearBuckling
from feastruct.solvers.naturalfrequency import NaturalFrequency
from feastruct.solvers.feasolve import SolverSettings

# ------------
# preprocessor
# ------------

# constants & lists
n = 20  # number of elements (make even to ensure point load in centre)
L = 10000  # length of arch
h = 3000  # height of arch
nodes = []  # list holding the node objects
elements = []  # list holding the element objects

# create 2d frame analysis object
analysis = FrameAnalysis2D()

# create materials and sections
steel = Steel()
section = Section(area=3230, ixx=23.6e6)

# create nodes
for i in range(n + 1):
    x = i * L / n
    y = h * np.sin(np.pi * x / L)
    nodes.append(analysis.create_node(coords=[x, y]))

# create beam elements
for i in range(n):
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
freedom_case.add_nodal_support(node=nodes[-1], val=0, dof=0)
freedom_case.add_nodal_support(node=nodes[-1], val=0, dof=1)

# add loads
load_case = cases.LoadCase()
load_case.add_nodal_load(node=nodes[int(n/2)], val=-1e3, dof=1)

# add analysis case
analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

# plot problem data
analysis.post.n_subdiv = 5  # reduce number of subdivisions for elements in post
analysis.post.plot_geom(analysis_case=analysis_case)

# -------------
# static solver
# -------------

settings = SolverSettings()
settings.linear_static.time_info = True

LinearStatic(analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings).solve()

# -----------
# static post
# -----------

analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=2500)
analysis.post.plot_frame_forces(analysis_case=analysis_case, axial=True)
analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)
analysis.post.plot_reactions(analysis_case=analysis_case)

# ---------------
# buckling solver
# ---------------

settings.linear_buckling.time_info = True

LinearBuckling(analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings).solve()

# -------------
# buckling post
# -------------

for i in range(settings.linear_buckling.num_modes):
    analysis.post.plot_buckling_results(analysis_case=analysis_case, buckling_mode=i)

# ----------------
# frequency solver
# ----------------

settings.natural_frequency.time_info = True

NaturalFrequency(
    analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings).solve()

# --------------
# frequency post
# --------------

for i in range(settings.natural_frequency.num_modes):
    analysis.post.plot_frequency_results(analysis_case=analysis_case, frequency_mode=i)
