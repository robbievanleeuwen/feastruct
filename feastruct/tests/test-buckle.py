import sys
sys.path.append('../')
from fea.frame2d import Frame2D
from solvers.linstatic import LinearStatic
from solvers.linbuckling import LinearBuckling
from post.post import PostProcessor

# N.B. using [N] and [mm]

analysis = Frame2D()

n_el = 12
L = 5000

# create nodes
for i in range(n_el + 1):
    analysis.add_node(id=i+1, coord=[0, i * L / n_el])

for i in range(n_el + 1, 2 * n_el + 1):
    analysis.add_node(id=i+1, coord=[(i - n_el) * L / n_el, L])

for i in range(2 * n_el + 1, 3 * n_el + 1):
    analysis.add_node(id=i+1, coord=[L, L - (i - 2 * n_el) * L / n_el])

# create elements
for i in range(3 * n_el):
    analysis.add_element(id=i+1, node_ids=[i+1, i+2], el_type='EB2', E=200e3,
                         A=1000, ixx=8.333e5)

fc1 = analysis.add_freedom_case(id=1)
fc1.add_nodal_support(node_id=1, val=0, dir=1)
fc1.add_nodal_support(node_id=1, val=0, dir=2)
fc1.add_nodal_support(node_id=3*n_el+1, val=0, dir=1)
fc1.add_nodal_support(node_id=3*n_el+1, val=0, dir=2)

lc1 = analysis.add_load_case(id=1)
lc1.add_nodal_load(node_id=n_el*1.5+1, val=-1e3, dir=2)

analysis.add_analysis_case(id=1, fc_id=1, lc_id=1)

post = PostProcessor(analysis, n_subdiv=5)
post.plot_geom(case_id=1)

LinearStatic(analysis, case_ids=[1]).solve()

post.plot_geom(case_id=1, deformed=True, def_scale=25)
post.plot_frame_forces(case_id=1, axial=True)
post.plot_frame_forces(case_id=1, shear=True)
post.plot_frame_forces(case_id=1, moment=True)

LinearBuckling(analysis, case_id=1).solve()

post.plot_buckling_eigenvector(case_id=1, buckling_mode=1)
post.plot_buckling_eigenvector(case_id=1, buckling_mode=2)
post.plot_buckling_eigenvector(case_id=1, buckling_mode=3)
post.plot_buckling_eigenvector(case_id=1, buckling_mode=4)
