import sys
sys.path.append('../')
from fea.frame2d import Frame2D
from solvers.linstatic import LinearStatic
from post.post import PostProcessor

# N.B. using [N] and [mm]

analysis = Frame2D()

analysis.add_node(id=1, coord=[0, 0])
analysis.add_node(id=2, coord=[1500, 0])
analysis.add_node(id=3, coord=[3000, 1500])

analysis.add_element(id=1, node_ids=[1, 2], el_type='EB2', E=200e3, A=3230,
                     ixx=23.6e6)
analysis.add_element(id=2, node_ids=[2, 3], el_type='EB2', E=200e3, A=3230,
                     ixx=23.6e6)

fc1 = analysis.add_freedom_case(id=1)
fc1.add_nodal_support(node_id=2, val=0, dir=1)
fc1.add_nodal_support(node_id=2, val=0, dir=2)
fc1.add_nodal_support(node_id=3, val=0, dir=2)

fc2 = analysis.add_freedom_case(id=2)
fc2.add_nodal_support(node_id=2, val=0, dir=1)
fc2.add_nodal_support(node_id=2, val=0, dir=2)
fc2.add_nodal_support(node_id=3, val=0, dir=2)
fc2.add_nodal_support(node_id=1, val=-10, dir=2)

fc3 = analysis.add_freedom_case(id=3)
fc3.add_nodal_support(node_id=2, val=0, dir=1)
fc3.add_nodal_support(node_id=2, val=0, dir=2)
fc3.add_nodal_support(node_id=3, val=0, dir=2)
fc3.add_nodal_support(node_id=1, val=-100, dir=2)
fc3.add_nodal_support(node_id=1, val=1, dir=3)

lc1 = analysis.add_load_case(id=1)
lc1.add_nodal_load(node_id=1, val=-1e4, dir=2)

lc2 = analysis.add_load_case(id=2)

analysis.add_analysis_case(id=1, fc_id=1, lc_id=1)
analysis.add_analysis_case(id=2, fc_id=2, lc_id=2)
analysis.add_analysis_case(id=3, fc_id=3, lc_id=2)

post = PostProcessor(analysis)
post.plot_geom(case_id=1)
post.plot_geom(case_id=2)
post.plot_geom(case_id=3)

LinearStatic(analysis, case_ids=[1, 2, 3]).solve()

post.plot_geom(case_id=1, deformed=True, def_scale=10)
post.plot_geom(case_id=2, deformed=True, def_scale=10)
post.plot_geom(case_id=3, deformed=True, def_scale=1)

post.plot_frame_forces(case_id=1, axial=True, shear=True)
post.plot_reactions(case_id=1)
