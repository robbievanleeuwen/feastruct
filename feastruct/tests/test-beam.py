from feastruct.fea.frame2d import Frame2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.post.post import PostProcessor

# N.B. using [N] and [mm]

analysis = Frame2D()

analysis.add_node(id=1, coord=[0, 0])
analysis.add_node(id=2, coord=[750, 0])
analysis.add_node(id=3, coord=[1500, 0])
analysis.add_node(id=4, coord=[2750, 0])
analysis.add_node(id=5, coord=[4000, 0])

analysis.add_element(id=1, node_ids=[1, 2], el_type='EB2', E=200e3, A=3230,
                     ixx=23.6e6)
analysis.add_element(id=2, node_ids=[2, 3], el_type='EB2', E=200e3, A=3230,
                     ixx=23.6e6)
analysis.add_element(id=3, node_ids=[3, 4], el_type='EB2', E=200e3, A=3230,
                     ixx=23.6e6)
analysis.add_element(id=4, node_ids=[4, 5], el_type='EB2', E=200e3, A=3230,
                     ixx=23.6e6)

fc1 = analysis.add_freedom_case(id=1)
fc1.add_nodal_support(node_id=1, val=0, dir=1)
fc1.add_nodal_support(node_id=1, val=0, dir=2)
fc1.add_nodal_support(node_id=3, val=0, dir=2)
fc1.add_nodal_support(node_id=5, val=0, dir=2)

lc1 = analysis.add_load_case(id=1)
lc1.add_nodal_load(node_id=2, val=-1e4, dir=2)
lc1.add_nodal_load(node_id=4, val=-1e4, dir=2)

analysis.add_analysis_case(id=1, fc_id=1, lc_id=1)

post = PostProcessor(analysis)
post.plot_geom(case_id=1)

LinearStatic(analysis, case_ids=[1]).solve()

post.plot_geom(case_id=1, deformed=True, def_scale=100)
post.plot_frame_forces(case_id=1, axial=True)
post.plot_frame_forces(case_id=1, shear=True)
post.plot_frame_forces(case_id=1, moment=True)
post.plot_reactions(case_id=1)
