# feastruct

[![Documentation Status](https://readthedocs.org/projects/feastruct/badge/?version=latest)](https://feastruct.readthedocs.io/en/latest/?badge=latest)

structural finite element analysis, the pythonic way.

*currently under development...*

```python
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.feasolve import SolverSettings

# ------------
# preprocessor
# ------------

# constants & lists
length = 5000  # length of the beam
udl = -10  # value of the udl

# everything starts with the analysis object
analysis = FrameAnalysis2D()

# materials and sections are objects
steel = Steel()
section = Section(area=3230, ixx=23.6e6)

# nodes are objects
node1 = analysis.create_node(coords=[0])
node2 = analysis.create_node(coords=[length])
node3 = analysis.create_node(coords=[2*length, 2500])

# and so are beams!
beam1 = analysis.create_element(
    el_type='EB2-2D', nodes=[node1, node2], material=steel, section=section
)
beam2 = analysis.create_element(
    el_type='EB2-2D', nodes=[node2, node3], material=steel, section=section
)

# boundary conditions are objects
freedom_case = cases.FreedomCase()
n1x = freedom_case.add_nodal_support(node=node1, val=0, dof=0)
n1y = freedom_case.add_nodal_support(node=node1, val=0, dof=1)
n2y = freedom_case.add_nodal_support(node=node2, val=0, dof=1)
n3y = freedom_case.add_nodal_support(node=node3, val=0, dof=1)

# so are loads!
load_case = cases.LoadCase()
udl1 = load_case.add_element_load(beam1.generate_udl(q=udl))
udl2 = load_case.add_element_load(beam2.generate_udl(q=udl))

# an analysis case relates a support case to a load case
analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

# ------
# solver
# ------

# you can easily change the solver settings
settings = SolverSettings()
settings.linear_static.time_info = True

# the linear static solver is an object and acts on the analysis object
LinearStatic(analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings).solve()

# ----
# post
# ----

# there are plenty of post processing options!
analysis.post.plot_geom(analysis_case=analysis_case)
analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=1e2)
analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True)
analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True)
analysis.post.plot_reactions(analysis_case=analysis_case)
```

## Current Capabilities:

### Pre-Processor
- [x] Python API
- [ ] Input File
- [ ] .dxf Import
- [ ] Triangular Meshing
- [ ] Structure Generator Functions

### Finite Element Analysis
- [x] 2D Frame
  - [x] Bar Element
  - [x] Euler Bernoulli Frame Element
  - [ ] Timoshenko Frame Element
- [ ] 3D Frame
  - [x] Bar Element
  - [ ] Euler Bernoulli Frame Element
  - [ ] Timoshenko Frame Element
- [ ] 2D Membrane (Plane Stress/Plane Strain)
  - [ ] 3-Noded Triangular Element
  - [ ] 6-Noded Triangular Element
- [ ] Plate Elements
- [ ] Shell Elements

### Element Formulations
- [x] Geometrically Linear
- [ ] Geometrically Non-Linear
- [ ] Material Non-Linear

### Loading/Restraints
- [x] Applied Loads
  - [x] Load Cases
  - [x] Nodal Loads
  - [x] Surface (Distributed) Loads
  - [ ] Body Loads
- [x] Restraints
  - [x] Freedom Cases
  - [x] Nodal Supports
  - [ ] Nodal Springs
  - [ ] Surface (Distributed) Supports
- [x] Analysis Cases

### Solvers
- [x] Linear Static Solver
- [ ] Non-Linear Static Solver
- [x] Linear Buckling Solver
- [x] Natural Frequency Solver
- [ ] Linear Dynamic Solver
- [ ] Non-Linear Dynamic Solver
- [x] Multi-Element Solvers

### Post-Processor
- [x] Structural Mesh and Boundary Conditions
- [x] Deformed Mesh
- [x] Reaction Forces
- [x] 2D Frame Actions (N, V, M)
- [x] Buckling Mode Shapes
- [x] Natural Frequency Mode Shapes
- [ ] Deformed Contour Plot
- [ ] Continuum Stress Contour Plot

### Additional Modules
- [ ] Optimisation
