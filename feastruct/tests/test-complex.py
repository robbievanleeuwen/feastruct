import sys
sys.path.append('../')
from fea.frame2d import Frame2D
from solvers.linstatic import LinearStatic
from post.post import PostProcessor

# N.B. using [N] and [mm]

analysis = Frame2D()

analysis.add_node(id=1, coord=[0, 0])
analysis.add_node(id=2, coord=[0, 4000])
analysis.add_node(id=3, coord=[3000, 0])
analysis.add_node(id=4, coord=[3000, 4000])
analysis.add_node(id=5, coord=[6000, 0])
analysis.add_node(id=6, coord=[6000, 4000])
analysis.add_node(id=7, coord=[9000, 2000])
analysis.add_node(id=8, coord=[9000, 6000])
analysis.add_node(id=9, coord=[10500, 6000])
analysis.add_node(id=10, coord=[12000, 2000])
analysis.add_node(id=11, coord=[15000, 5000])

analysis.add_element(id=1, node_ids=[1, 2])
analysis.add_element(id=2, node_ids=[2, 4])
analysis.add_element(id=3, node_ids=[3, 4])
analysis.add_element(id=4, node_ids=[4, 6])
analysis.add_element(id=5, node_ids=[5, 6])
analysis.add_element(id=6, node_ids=[6, 7])
analysis.add_element(id=7, node_ids=[6, 8])
analysis.add_element(id=8, node_ids=[7, 10])
analysis.add_element(id=9, node_ids=[8, 9])
analysis.add_element(id=10, node_ids=[9, 10])
analysis.add_element(id=11, node_ids=[10, 11])
analysis.add_element(id=12, node_ids=[9, 11])

fc1 = analysis.add_freedom_case(id=1)
fc1.add_nodal_support(node_id=1, val=0, dir=1)
fc1.add_nodal_support(node_id=1, val=0, dir=2)
fc1.add_nodal_support(node_id=1, val=0, dir=3)
fc1.add_nodal_support(node_id=2, val=0, dir=1)
fc1.add_nodal_support(node_id=3, val=0, dir=1)
fc1.add_nodal_support(node_id=3, val=0, dir=2)
fc1.add_nodal_support(node_id=4, val=0, dir=3)
fc1.add_nodal_support(node_id=5, val=0, dir=2)
fc1.add_nodal_support(node_id=7, val=0, dir=2)
fc1.add_nodal_support(node_id=7, val=0, dir=3)
fc1.add_nodal_support(node_id=8, val=0, dir=1)
fc1.add_nodal_support(node_id=8, val=0, dir=3)
fc1.add_nodal_support(node_id=11, val=0, dir=1)
fc1.add_nodal_support(node_id=11, val=0, dir=2)
fc1.add_nodal_support(node_id=6, val=0.1, dir=1)
fc1.add_nodal_support(node_id=9, val=-10, dir=2)
fc1.add_nodal_support(node_id=5, val=-1, dir=3)

lc1 = analysis.add_load_case(id=1)
lc1.add_nodal_load(node_id=2, val=-1, dir=2)
lc1.add_nodal_load(node_id=10, val=1, dir=1)
lc1.add_nodal_load(node_id=11, val=-1, dir=3)

analysis.add_analysis_case(id=1, fc_id=1, lc_id=1)

post = PostProcessor(analysis)
post.plot_geom(case_id=1)

solver = LinearStatic(analysis, case_ids=[1]).solve()

post.plot_geom(case_id=1, deformed=True, def_scale=0.005)
post.plot_reactions(case_id=1)
