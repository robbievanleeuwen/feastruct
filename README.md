# feastruct
a python package for structural finite element analysis.

*currently under development...*

## Current Capabilities:

### Pre-Processor
- [x] Python Interface
- [ ] Input File
- [ ] .dxf Import
- [ ] Triangular Meshing
- [ ] Structure Generator Functions

### Finite Element Analysis
- [x] 2D Frame
  - [x] Euler Bernoulli Frame Element
  - [ ] Timoshenko Frame Element
  - [ ] Truss Element
- [ ] 2D Membrane (Plane Stress/Plane Strain)
  - [ ] 3-Noded Triangular Element
  - [ ] 6-Noded Triangular Element
- [ ] Plate Elements
- [ ] Shell Elements
- [ ] 3D Frame Elements

### Element Formulations
- [x] Geometrically Linear
- [ ] Geometrically Non-Linear
- [ ] Material Non-Linearity

### Loading/Restraints
- [x] Applied Loads
  - [x] Nodal Loads
  - [ ] Surface (Distributed) Loads
  - [ ] Body Loads
  - [x] Load Cases
- [x] Restraints
  - [x] Nodal Supports
  - [ ] Surface Supports
  - [ ] Nodal Springs
  - [x] Freedom Cases
- [x] Analysis Cases

### Solvers
- [x] Linear Static Solver
- [ ] Non-Linear Static Solver
- [x] Linear Buckling Solver
- [x] Natural Frequency Solver
- [ ] Linear Dynamic Solver
- [ ] Non-Linear Dynamic Solver

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
