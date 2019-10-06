import time
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
from feastruct.pre.material import Steel
from feastruct.pre.section import Section
import feastruct.fea.cases as cases
from feastruct.fea.frame_analysis import FrameAnalysis2D
from feastruct.solvers.linstatic import LinearStatic
from feastruct.solvers.naturalfrequency import NaturalFrequency
from feastruct.solvers.feasolve import SolverSettings


def preprocessor(length, depth, panels, freq):
    w = 10  # load per unit length [N/mm]
    dx = length / panels  # length of a truss panel

    # create 2d frame analysis object
    analysis = FrameAnalysis2D()

    # create materials and sections
    steel = Steel()
    section_deck = Section(area=3230, ixx=23.6e6)  # 200UB25
    section_vertical = Section(area=3000, ixx=3.91e6)  # 100x9 SHS
    section_diagonal = Section(area=314)  # 20 dia rod
    section_tie = Section(area=1018)  # 36 dia rod

    # create deck nodes
    nodes_deck = []

    for i in range(panels + 1):
        nodes_deck.append(analysis.create_node(coords=[i*dx]))

    # create tie nodes - first node coincides with first deck node
    nodes_tie = [nodes_deck[0]]

    for i in range(panels - 1):
        # fit parabola
        x = (i + 1) * dx
        y = 4 * depth / length / length * x * x - 4 * depth / length * x

        nodes_tie.append(analysis.create_node(coords=[x, y]))

    # last node coincides with last deck node
    nodes_tie.append(nodes_deck[-1])

    # create deck (beam) elements
    elements_deck = []

    for i in range(panels):
        elements_deck.append(analysis.create_element(
            el_type='EB2-2D',
            nodes=[nodes_deck[i], nodes_deck[i+1]],
            material=steel,
            section=section_deck
        ))

    # create tie (bar) elements
    elements_tie = []

    for i in range(panels):
        elements_tie.append(analysis.create_element(
            el_type='Bar2-2D',
            nodes=[nodes_tie[i], nodes_tie[i+1]],
            material=steel,
            section=section_tie
        ))

    # create vertical (beam) elements
    elements_vertical = []

    for i in range(panels - 1):
        elements_vertical.append(analysis.create_element(
            el_type='EB2-2D',
            nodes=[nodes_deck[i+1], nodes_tie[i+1]],
            material=steel,
            section=section_vertical
        ))

    # create diagonal (bar) elements
    elements_diaongal = []

    for i in range(panels - 2):
        elements_diaongal.append(analysis.create_element(
            el_type='Bar2-2D',
            nodes=[nodes_deck[i+1], nodes_tie[i+2]],
            material=steel,
            section=section_diagonal
        ))

        elements_diaongal.append(analysis.create_element(
            el_type='Bar2-2D',
            nodes=[nodes_deck[i+2], nodes_tie[i+1]],
            material=steel,
            section=section_diagonal
        ))

    # add supports
    freedom_case = cases.FreedomCase()
    freedom_case.add_nodal_support(node=nodes_deck[0], val=0, dof=0)
    freedom_case.add_nodal_support(node=nodes_deck[0], val=0, dof=1)
    freedom_case.add_nodal_support(node=nodes_deck[-1], val=0, dof=1)

    # add loads
    load_case = cases.LoadCase()

    # if frequency analysis, don't bother adding load
    if not freq:
        for el in elements_deck:
            load_case.add_element_load(el.generate_udl(q=-w))

    # add analysis case
    analysis_case = cases.AnalysisCase(freedom_case=freedom_case, load_case=load_case)

    return (analysis, analysis_case)


# truss geometry
length = 20000  # length of the truss [mm]
depth = 2500  # depth of the truss at midspan [mm]
panels = 6  # number of truss panels

(analysis, analysis_case) = preprocessor(length, depth, panels, False)

# -------------
# static solver
# -------------

settings = SolverSettings()
settings.linear_static.time_info = True

LinearStatic(analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings).solve()

# -----------
# static post
# -----------

analysis.post.plot_geom(analysis_case=analysis_case, deformed=True, def_scale=25)
analysis.post.plot_frame_forces(analysis_case=analysis_case, axial=True, scale=0.05)
analysis.post.plot_frame_forces(analysis_case=analysis_case, shear=True, scale=0.05)
analysis.post.plot_frame_forces(analysis_case=analysis_case, moment=True, scale=0.05)

# ----------------
# frequency solver
# ----------------

settings.natural_frequency.time_info = True

solver = NaturalFrequency(
    analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings)
solver.solve()

# --------------
# frequency post
# --------------

for i in range(settings.natural_frequency.num_modes):
    analysis.post.plot_frequency_results(analysis_case=analysis_case, frequency_mode=i)

# ------------
# optimisation
# ------------

# for varying spans, find required depths such that first mode equals 8 Hz


def frequency_analysis(depth, length, opt):
    """objective function"""

    panels = 6  # number of truss panels

    # build preprocessor
    (analysis, analysis_case) = preprocessor(length, depth, panels, True)

    # solve only for the first mode
    settings = SolverSettings()
    settings.natural_frequency.num_modes = 1

    solver = NaturalFrequency(
        analysis=analysis, analysis_cases=[analysis_case], solver_settings=settings)
    solver.solve()

    # extract first frequency from any element within the analysis
    (f, v) = analysis.elements[0].get_frequency_results(analysis_case=analysis_case)

    # return the difference between the frequency and 8 Hz
    if opt:
        return f - 8
    else:
        return (analysis, analysis_case)


# vary spans between 10 m and 40 m
lengths = np.linspace(10e3, 40e3, 11)
depths = []
print("")

for length in lengths:
    # search between depths of 50 mm and 8000 mm
    # N.B. be careful with initial search bounds (may introduce lower order local modes)
    start_time = time.time()
    (d, r) = optimize.brentq(
        f=frequency_analysis, a=50, b=8000, args=(length, True), full_output=True)
    sol_time = time.time() - start_time

    print("Optimisation converged in {0:.5f} seconds after {1} iterations.".format(
        sol_time, r.iterations))
    print("Depth of {0:.2f} mm gives a natural frequency of 8 Hz for a span of {1:.2f} mm.".format(
        d, length))

    depths.append(d)

fig = plt.figure()
plt.plot(lengths, depths, 'kx-')
plt.title('Bowstring Truss Depth for $f_1$ = 8 Hz')
plt.xlabel('Truss Span [mm]')
plt.ylabel('Truss Depth [mm]')
plt.show()

(analysis, analysis_case) = frequency_analysis(depth=depths[-1], length=lengths[-1], opt=False)
analysis.post.plot_frequency_results(analysis_case)
