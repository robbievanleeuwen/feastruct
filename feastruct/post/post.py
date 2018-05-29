import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from fea.exceptions import FEAInputError


class PostProcessor:
    """docstring for PostProcessor.
    """

    def __init__(self, analysis):
        """asldkjasld
        """

        self.analysis = analysis
        self.n_subdiv = 50

    def plot_geom(self, case_id, undeformed=True, deformed=False, def_scale=1):
        """askldjasld

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        # get analysis case
        try:
            analysis_case = self.analysis.find_analysis_case(case_id)
        except FEAInputError as error:
            print(error)
            sys.exit(1)

        (fig, ax) = plt.subplots()

        for el in self.analysis.elements:
            if deformed:
                el.plot_deformed_element(ax, case_id, self.n_subdiv, def_scale)
                if undeformed:
                    el.plot_element(ax, linestyle='--', linewidth=1, marker='')
            else:
                el.plot_element(ax)

        # set initial plot limits
        (xmin, xmax, ymin, ymax) = self.analysis.get_node_lims()
        ax.set_xlim(xmin-1e-12, xmax)
        ax.set_ylim(ymin-1e-12, ymax)

        # get 2% of the maxmimum dimension
        small = 0.02*max(xmax-xmin, ymax-ymin)

        # generate lists of nodal supports and imposed displacements
        support_node_list = []
        imposed_disp_list = []
        max_disp = 0

        for support in analysis_case.freedom_case.items:
            # check that there is no imposed displacement
            if support.val == 0:
                support.node.fixity[support.dir-1] = 1

                if support.node not in support_node_list:
                    support_node_list.append(support)
            # if there is an imposed displacement
            else:
                imposed_disp_list.append(support)
                if support.dir in (1, 2):
                    max_disp = max(max_disp, abs(support.val))

        # plot supports
        for support in support_node_list:
            support.plot_support(ax, max_disp, small, self.get_support_angle,
                                 case_id, deformed, def_scale)

        # plot imposed displacements
        for imposed_disp in imposed_disp_list:
            if imposed_disp.dir in (1, 2):
                imposed_disp.plot_imposed_disp(
                    ax, max_disp, small, self.get_support_angle, case_id,
                    deformed, def_scale)
            elif imposed_disp.dir == 3:
                imposed_disp.plot_imposed_rot(
                    ax, small, self.get_support_angle, case_id, deformed,
                    def_scale)

        # find max force
        max_force = 0

        for load in analysis_case.load_case.items:
            if load.dir == 1 or load.dir == 2:
                max_force = max(max_force, abs(load.val))

        # plot loads
        for load in analysis_case.load_case.items:
            load.plot_load(ax, max_force, small, self.get_support_angle,
                           case_id, deformed, def_scale)

        # plot layout
        plt.axis('tight')
        ax.set_xlim(wide_lim(ax.get_xlim()))
        ax.set_ylim(wide_lim(ax.get_ylim()))

        limratio = np.diff(ax.get_ylim())/np.diff(ax.get_xlim())

        if limratio < 0.5:
            ymid = np.mean(ax.get_ylim())
            ax.set_ylim(ymid + (ax.get_ylim() - ymid) * 0.5 / limratio)
        elif limratio > 1:
            xmid = np.mean(ax.get_xlim())
            ax.set_xlim(xmid + (ax.get_xlim() - xmid) * limratio)

        ax.set_aspect(1)
        plt.box(on=None)
        plt.show()

        return ax

    def get_support_angle(self, node, prefer_dir=None):
        """alskdjaklsd

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        # find angles to connected elements
        phi = []
        num_el = 0

        # loop through each element in the mesh
        for el in self.analysis.elements:
            # if the current element is connected to the node
            if node in el.nodes:
                num_el += 1
                # loop through all the nodes connected to the element
                for el_node in el.nodes:
                    # if the node is not the node in question
                    if el_node is not node:
                        dx = [el_node.x - node.x, el_node.y - node.y]
                        phi.append(np.arctan2(dx[1], dx[0]) / np.pi * 180)

        phi.sort()
        phi.append(phi[0] + 360)
        i0 = np.argmax(np.diff(phi))
        angle = (phi[i0] + phi[i0+1]) / 2 + 180

        if prefer_dir is not None:
            if prefer_dir == 2:
                if max(np.sin([phi[i0] * np.pi / 180,
                               phi[i0+1] * np.pi / 180])) > -0.1:
                    angle = 90
            elif prefer_dir == 1:
                if max(np.cos([phi[i0] * np.pi / 180,
                               phi[i0+1] * np.pi / 180])) > -0.1:
                    angle = 0

        return (angle, num_el)


def wide_lim(x):
    x2 = max(x)
    x1 = min(x)
    dx = x2-x1
    f = 0.02
    return (x1-f*dx, x2+f*dx)
