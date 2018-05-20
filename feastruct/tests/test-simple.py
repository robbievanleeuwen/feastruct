import sys
sys.path.append('../')
from fea.frame2d import Frame2D
from solvers.linstatic import LinearStatic
from post.post import PostProcessor

# N.B. using [N] and [mm]

analysis = Frame2D(analysis_type='linear')

analysis.add_node(id=1, coord=[0, 0])
analysis.add_node(id=2, coord=[1500, 0])
analysis.add_node(id=3, coord=[3000, 1500])

analysis.add_element(id=1, node_ids=[1, 2], el_type='TB2', E=200e3, A=3230,
                     ixx=23.6e6, G=1, A_s=1e6)
analysis.add_element(id=2, node_ids=[2, 3], el_type='TB2', E=200e3, A=3230,
                     ixx=23.6e6, G=1, A_s=1e6)

analysis.add_support(node_id=2, val=0, dir=1)
analysis.add_support(node_id=2, val=0, dir=2)
analysis.add_support(node_id=3, val=0, dir=2)

analysis.add_nodal_load(node_id=1, val=-1e4, dir=2)

post = PostProcessor(analysis)
post.plot_geom()

solver = LinearStatic(analysis).solve()
