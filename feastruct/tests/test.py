import sys
sys.path.append('../')
from fea.frame2d import Frame2D
from solvers.linstatic import LinearStatic

analysis = Frame2D(analysis_type='linear')

analysis.add_node(id=1, coord=[0, 0])
analysis.add_node(id=2, coord=[1, 0])
analysis.add_node(id=3, coord=[3, 5])

analysis.add_element(id=1, node_ids=[1, 2], el_type='EB2')
analysis.add_element(id=2, node_ids=[2, 3], el_type='TB2')

analysis.add_support(node_id=1, val=0, dir=1)
analysis.add_support(node_id=1, val=0, dir=2)
analysis.add_support(node_id=3, val=0, dir=2)

analysis.add_nodal_load(node_id=2, val=-1, dir=2)

solver = LinearStatic(analysis)
solver.solve()
