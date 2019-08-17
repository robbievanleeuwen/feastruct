# feastruct
a python package for structural finite element analysis.

*currently under development...*

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
