import numpy as np
import matplotlib.pyplot as plt
from feastruct.fea.exceptions import FEAInputError


class PostProcessor:
    """Class for post processing methods.

    This class provides some post-processing methods for a particular analysis
    that can be used to visualise tthe structural geometry and the finite
    element analysis results.

    :cvar analysis: Analysis object for post-processing
    :vartype analysis: :class:`feastruct.fea.fea.fea`
    :cvar int n_subdiv: Number of subdivisions used to discretise frame
        elements in post-processing, such that higher order shape functions can
        be realised
    """

    def __init__(self, analysis, n_subdiv=20):
        """Inits the fea class.

        :param analysis: Analysis object for post-processing
        :type analysis: :class:`feastruct.fea.fea.fea`
        :param int n_subdiv: Number of subdivisions used to discretise frame
            elements in post-processing
        """

        self.analysis = analysis
        self.n_subdiv = n_subdiv

    def plot_geom(self, case_id, ax=None, supports=True, loads=True,
                  undeformed=True, deformed=False, def_scale=1, dashed=False):
        """Method used to plot the structural mesh in the undeformed and/or
        deformed state. If no axes object is provided, a new axes object is
        created. N.B. this method is adopted from the MATLAB code by F.P. van
        der Meer: plotGeom.m.

        :param int case_id: Unique case id
        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param bool supports: Whether or not the freedom case supports are
            rendered
        :param bool loads: Whether or not the load case loads are rendered
        :param bool undeformed: Whether or not the undeformed structure is
            plotted
        :param bool deformed: Whether or not the deformed structure is plotted
        :param float def_scale: Deformation scale used for plotting the
            deformed structure
        :param bool dashed: Whether or not to plot the structure with dashed
            lines if only the undeformed structure is to be plotted
        """

        # get analysis case
        try:
            analysis_case = self.analysis.find_analysis_case(case_id)
        except FEAInputError as error:
            print(error)

        if ax is None:
            (fig, ax) = plt.subplots()

        for el in self.analysis.elements:
            if deformed:
                el.plot_deformed_element(ax, case_id, self.n_subdiv, def_scale)
                if undeformed:
                    el.plot_element(ax, linestyle='--', linewidth=1, marker='')
            else:
                if dashed:
                    el.plot_element(ax, linestyle='--', linewidth=1, marker='')
                else:
                    el.plot_element(ax)

        # set initial plot limits
        (xmin, xmax, ymin, ymax) = self.analysis.get_node_lims()
        ax.set_xlim(xmin-1e-12, xmax)
        ax.set_ylim(ymin-1e-12, ymax)

        # get 2% of the maxmimum dimension
        small = 0.02 * max(xmax-xmin, ymax-ymin)

        if supports:
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

        if supports:
            # plot supports
            for support in support_node_list:
                support.plot_support(ax, small, self.get_support_angle,
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

        if loads:
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
        ax.set_xlim(self.wide_lim(ax.get_xlim()))
        ax.set_ylim(self.wide_lim(ax.get_ylim()))

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

    def plot_reactions(self, case_id):
        """Method used to generate a plot of the reaction forces.

        :param int case_id: Unique case id
        """

        (fig, ax) = plt.subplots()

        # get analysis case
        try:
            analysis_case = self.analysis.find_analysis_case(case_id)
        except FEAInputError as error:
            print(error)

        # get size of structure
        (xmin, xmax, ymin, ymax) = self.analysis.get_node_lims()

        # determine maximum reaction force
        max_reaction = 0

        for support in analysis_case.freedom_case.items:
            if support.dir in (1, 2):
                try:
                    reaction = support.get_reaction(case_id)
                except FEAInputError as error:
                    print(error)
                max_reaction = max(max_reaction, abs(reaction))

        small = 0.02 * max(xmax-xmin, ymax-ymin)

        # plot reactions
        for support in analysis_case.freedom_case.items:
            support.plot_reaction(ax, max_reaction, small,
                                  self.get_support_angle, case_id)

        # plot the undeformed structure
        self.plot_geom(case_id, ax=ax, supports=False)

    def plot_frame_forces(self, case_id, axial=False, shear=False,
                          moment=False, scale=0.1):
        """Method used to plot internal frame actions resulting from the
        analysis case case_id.

        :param int case_id: Unique case id
        :param bool axial: Whether or not the axial force diagram is displayed
        :param bool shear: Whether or not the shear force diagram is displayed
        :param bool moment: Whether or not the bending moment diagram is
            displayed
        :param float scale: Scale used for plotting internal force diagrams.
            Corresponds to the fraction of the window that the largest action
            takes up
        """

        # TODO: check that analysis object is Frame2D

        (fig, ax) = plt.subplots()

        # get size of structure
        (xmin, xmax, ymin, ymax) = self.analysis.get_node_lims()

        # determine maximum forces
        max_axial = 0
        max_shear = 0
        max_moment = 0

        # loop throuh each element to get max forces
        for el in self.analysis.elements:
            try:
                f_int = el.get_fint(case_id)
                if axial:
                    max_axial = max(max_axial, abs(f_int.N1), abs(f_int.N2))
                if shear:
                    max_shear = max(max_shear, abs(f_int.V1), abs(f_int.V2))
                if moment:
                    max_moment = max(max_moment, abs(f_int.M1), abs(f_int.M2))
            except FEAInputError as error:
                print(error)

        scale_axial = scale * max(xmax - xmin, ymax - ymin) / max(
            max_axial, 1e-8)
        scale_shear = scale * max(xmax - xmin, ymax - ymin) / max(
            max_shear, 1e-8)
        scale_moment = scale * max(xmax - xmin, ymax - ymin) / max(
            max_moment, 1e-8)

        # loop throgh each element to plot the forces
        for el in self.analysis.elements:
            if axial:
                el.plot_axial_force(ax, case_id, scale_axial)
            if shear:
                el.plot_shear_force(ax, case_id, scale_shear)
            if moment:
                el.plot_bending_moment(ax, case_id, scale_moment)

        # plot the undeformed structure
        self.plot_geom(case_id, ax=ax)

    def plot_buckling_eigenvector(self, case_id, buckling_mode=1):
        """Method used to plot a buckling eigenvector. The undeformed structure
        is plotted with a dahsed line.

        :param int case_id: Unique case id
        :param int buckling_mode: Buckling mode to plot
        """

        (fig, ax) = plt.subplots()

        # set initial plot limits
        (xmin, xmax, ymin, ymax) = self.analysis.get_node_lims()

        # determine max eigenvector displacement value (ignore rotation)
        max_v = 0

        for el in self.analysis.elements:
            (w, v_el) = el.get_buckling_results(case_id, buckling_mode)
            max_v = max(max_v, abs(v_el[0, 0]), abs(v_el[0, 1]),
                        abs(v_el[1, 0]), abs(v_el[1, 1]))

        # determine plot scale
        scale = 0.1 * max(xmax - xmin, ymax - ymin) / max_v

        # plot eigenvectors
        for el in self.analysis.elements:
            (_, v_el) = el.get_buckling_results(case_id, buckling_mode)
            el.plot_deformed_element(ax, case_id, self.n_subdiv, scale, v_el)

        # plot the load factor (eigenvalue)
        ax.set_title("Load Factor for Mode {:d}: {:.4e}".format(
            buckling_mode, w), size=10)

        # plot the undeformed structure
        self.plot_geom(case_id, ax=ax, dashed=True)

    def plot_frequency_eigenvector(self, case_id, frequency_mode=1):
        """Method used to plot a natural frequency eigenvector. The undeformed
        structure is plotted with a dahsed line.

        :param int case_id: Unique case id
        :param int frequency_mode: Frequency mode to plot
        """

        (fig, ax) = plt.subplots()

        # set initial plot limits
        (xmin, xmax, ymin, ymax) = self.analysis.get_node_lims()

        # determine max eigenvector displacement value (ignore rotation)
        max_v = 0

        for el in self.analysis.elements:
            (w, v_el) = el.get_frequency_results(case_id, frequency_mode)
            max_v = max(max_v, abs(v_el[0, 0]), abs(v_el[0, 1]),
                        abs(v_el[1, 0]), abs(v_el[1, 1]))

        # determine plot scale
        scale = 0.1 * max(xmax - xmin, ymax - ymin) / max_v

        # plot eigenvectors
        for el in self.analysis.elements:
            (_, v_el) = el.get_frequency_results(case_id, frequency_mode)
            el.plot_deformed_element(ax, case_id, self.n_subdiv, scale, v_el)

        # plot the natural frequency
        ax.set_title("Natural Frequency for Mode {:d}: {:.4e} Hz".format(
            frequency_mode, w), size=10)

        # plot the undeformed structure
        self.plot_geom(case_id, ax=ax, dashed=True)

    def get_support_angle(self, node, prefer_dir=None):
        """Given a node object, returns the optimal angle to plot a support.
        Essentially finds the average angle of the connected elements and
        considers a preferred plotting direction. N.B. this method is adopted
        from the MATLAB code by F.P. van der Meer: plotGeom.m.

        :param node: Node object
        :type node: :class:`feastruct.fea.node.node`
        :param int prefer_dir: Preferred direction to plot the support, where 1
            corresponds to the x-axis and 2 corresponds to the y-axis
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

    def wide_lim(self, x):
        """Returns a tuple corresponding to the axis limits (x) stretched by 2%
        on either side.

        :param x: List containing axis limits e.g. [xmin, xmax]
        :type x: list[float, float]
        :return: Stretched axis limits (x1, x2)
        :rtype: tuple(float, float)
        """

        x2 = max(x)
        x1 = min(x)
        dx = x2-x1
        f = 0.02
        return (x1-f*dx, x2+f*dx)
