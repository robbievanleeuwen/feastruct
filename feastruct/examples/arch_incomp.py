import numpy as np
from feastruct.fea.frame2d import Frame2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.linbuckling import LinearBuckling
from feastruct.solvers.naturalfrequency import NaturalFrequency
from feastruct.post.post import PostProcessor

# N.B. using [N] and [mm]

analysis = Frame2D()

n = 20  # number of elements
L = 10000  # length of arch
h = 3000  # height of arch

for i in range(n + 1):
    x = i * L / n
    y = h * np.sin(np.pi * x / L)
    analysis.add_node(id=i+1, coord=[x, y])

for i in range(n):
    analysis.add_element(id=i+1, node_ids=[i+1, i+2], el_type='EB2', E=200e3,
                         A=3230, ixx=23.6e6, rho=7.85e-9)

fc = analysis.add_freedom_case(id=1)
fc.add_nodal_support(node_id=1, val=0, dir=1)
fc.add_nodal_support(node_id=1, val=0, dir=2)
fc.add_nodal_support(node_id=n+1, val=0, dir=1)
fc.add_nodal_support(node_id=n+1, val=0, dir=2)

lc = analysis.add_load_case(id=1)
lc.add_nodal_load(node_id=n/2+1, val=-1e3, dir=2)

analysis.add_analysis_case(id=1, fc_id=1, lc_id=1)

post = PostProcessor(analysis, n_subdiv=5)
post.plot_geom(case_id=1)

LinearStatic(analysis, case_ids=[1]).solve()

post.plot_geom(case_id=1, deformed=True, def_scale=2500)

post.plot_frame_forces(case_id=1, axial=True)
post.plot_frame_forces(case_id=1, shear=True)
post.plot_frame_forces(case_id=1, moment=True)
post.plot_reactions(case_id=1)

LinearBuckling(analysis, case_id=1).solve()

post.plot_buckling_eigenvector(case_id=1, buckling_mode=1)
post.plot_buckling_eigenvector(case_id=1, buckling_mode=2)
post.plot_buckling_eigenvector(case_id=1, buckling_mode=3)
post.plot_buckling_eigenvector(case_id=1, buckling_mode=4)

NaturalFrequency(analysis, case_id=1).solve()

post.plot_frequency_eigenvector(case_id=1, frequency_mode=1)
post.plot_frequency_eigenvector(case_id=1, frequency_mode=2)
post.plot_frequency_eigenvector(case_id=1, frequency_mode=3)
post.plot_frequency_eigenvector(case_id=1, frequency_mode=4)
