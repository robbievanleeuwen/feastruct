import numpy as np
import matplotlib.pyplot as plt


class PostProcessor:
    """Class for post processing methods.

    This class provides some post-processing methods for a particular analysis that can be used to
    visualise tthe structural geometry and the finite element analysis results.

    :cvar analysis: Analysis object for post-processing
    :vartype analysis: :class:`~feastruct.fea.fea.fea`
    :cvar int n_subdiv: Number of subdivisions used to discretise frame
        elements in post-processing, such that higher order shape functions can
        be realised
    """

    def __init__(self, analysis, n_subdiv=20):
        """Inits the fea class.

        :param analysis: Analysis object for post-processing
        :type analysis: :class:`~feastruct.fea.fea.FiniteElementAnalysis`
        :param int n_subdiv: Number of subdivisions used to discretise frame elements in
            post-processing
        """

        self.analysis = analysis
        self.n_subdiv = n_subdiv

    def plot_geom(self, analysis_case, ax=None, supports=True, loads=True, undeformed=True,
                  deformed=False, def_scale=1, dashed=False):
        """Method used to plot the structural mesh in the undeformed and/or deformed state. If no
        axes object is provided, a new axes object is created. N.B. this method is adopted from the
        MATLAB code by F.P. van der Meer: plotGeom.m.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param bool supports: Whether or not the freedom case supports are rendered
        :param bool loads: Whether or not the load case loads are rendered
        :param bool undeformed: Whether or not the undeformed structure is plotted
        :param bool deformed: Whether or not the deformed structure is plotted
        :param float def_scale: Deformation scale used for plotting the deformed structure
        :param bool dashed: Whether or not to plot the structure with dashed lines if only the
            undeformed structure is to be plotted
        """

        if ax is None:
            (fig, ax) = plt.subplots()

        for el in self.analysis.elements:
            if deformed:
                el.plot_deformed_element(
                    ax=ax, analysis_case=analysis_case, n=self.n_subdiv, def_scale=def_scale)
                if undeformed:
                    el.plot_element(ax=ax, linestyle='--', linewidth=1, marker='')
            else:
                if dashed:
                    el.plot_element(ax=ax, linestyle='--', linewidth=1, marker='')
                else:
                    el.plot_element(ax=ax)

        # set initial plot limits
        (xmin, xmax, ymin, ymax, _, _) = self.analysis.get_node_lims()
        ax.set_xlim(xmin-1e-12, xmax)
        ax.set_ylim(ymin-1e-12, ymax)

        # get 2% of the maxmimum dimension
        small = 0.02 * max(xmax-xmin, ymax-ymin)

        if supports:
            # generate list of supports and imposed displacements for unique nodes
            node_supports = []
            node_imposed_disps = []
            max_disp = 0

            # loop through supports
            for support in analysis_case.freedom_case.items:
                # if there is no imposed displacement, the support is fixed
                if support.val == 0:
                    # check to see if the node hasn't already been added
                    if support.node not in node_supports:
                        node_supports.append(support)

                # if there is an imposed displacement
                else:
                    node_imposed_disps.append(support)
                    if support.dof in [0, 1]:
                        max_disp = max(max_disp, abs(support.val))

        if supports:
            # plot supports
            for support in node_supports:
                support.plot_support(
                    ax=ax, small=small, get_support_angle=self.get_support_angle,
                    analysis_case=analysis_case, deformed=deformed, def_scale=def_scale)

            # plot imposed displacements
            for imposed_disp in node_imposed_disps:
                if imposed_disp.dof in [0, 1]:
                    imposed_disp.plot_imposed_disp(
                        ax=ax, max_disp=max_disp, small=small,
                        get_support_angle=self.get_support_angle, analysis_case=analysis_case,
                        deformed=deformed, def_scale=def_scale)
                elif imposed_disp.dof == 5:
                    imposed_disp.plot_imposed_rot(
                        ax=ax, small=small, get_support_angle=self.get_support_angle,
                        analysis_case=analysis_case, deformed=deformed, def_scale=def_scale)

        if loads:
            # find max force
            max_force = 0

            for load in analysis_case.load_case.items:
                if load.dof in [0, 1]:
                    max_force = max(max_force, abs(load.val))

            # plot loads
            for load in analysis_case.load_case.items:
                load.plot_load(
                    ax=ax, max_force=max_force, small=small,
                    get_support_angle=self.get_support_angle, analysis_case=analysis_case,
                    deformed=deformed, def_scale=def_scale)

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

    def plot_reactions(self, analysis_case):
        """Method used to generate a plot of the reaction forces.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        (fig, ax) = plt.subplots()

        # get size of structure
        (xmin, xmax, ymin, ymax, _, _) = self.analysis.get_node_lims()

        # determine maximum reaction force
        max_reaction = 0

        for support in analysis_case.freedom_case.items:
            if support.dof in [0, 1]:
                reaction = support.get_reaction(analysis_case=analysis_case)
                max_reaction = max(max_reaction, abs(reaction))

        small = 0.02 * max(xmax-xmin, ymax-ymin)

        # plot reactions
        for support in analysis_case.freedom_case.items:
            support.plot_reaction(
                ax=ax, max_reaction=max_reaction, small=small,
                get_support_angle=self.get_support_angle, analysis_case=analysis_case)

        # plot the undeformed structure
        self.plot_geom(analysis_case=analysis_case, ax=ax, supports=False)

    def plot_frame_forces(self, analysis_case, axial=False, shear=False, moment=False, scale=0.1):
        """Method used to plot internal frame actions resulting from the analysis case.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param bool axial: Whether or not the axial force diagram is displayed
        :param bool shear: Whether or not the shear force diagram is displayed
        :param bool moment: Whether or not the bending moment diagram is displayed
        :param float scale: Scale used for plotting internal force diagrams. Corresponds to the
            fraction of the window that the largest action takes up
        """

        # TODO: implement for arbitrary number of stations - label ends and min, max

        (fig, ax) = plt.subplots()

        # get size of structure
        (xmin, xmax, ymin, ymax, _, _) = self.analysis.get_node_lims()

        # determine maximum forces
        max_axial = 0
        max_shear = 0
        max_moment = 0

        # loop throuh each element to get max forces
        for el in self.analysis.elements:
            if axial:
                afd = el.get_afd(n=2, analysis_case=analysis_case)
                max_axial = max(max_axial, abs(afd[0]), abs(afd[1]))
            if shear:
                sfd = el.get_sfd(n=2, analysis_case=analysis_case)
                max_shear = max(max_shear, abs(sfd[0]), abs(sfd[1]))
            if moment:
                bmd = el.get_bmd(n=2, analysis_case=analysis_case)
                max_moment = max(max_moment, abs(bmd[0]), abs(bmd[1]))

        scale_axial = scale * max(xmax - xmin, ymax - ymin) / max(max_axial, 1e-8)
        scale_shear = scale * max(xmax - xmin, ymax - ymin) / max(max_shear, 1e-8)
        scale_moment = scale * max(xmax - xmin, ymax - ymin) / max(max_moment, 1e-8)

        # loop throgh each element to plot the forces
        for el in self.analysis.elements:
            if axial:
                el.plot_axial_force(ax=ax, analysis_case=analysis_case, scalef=scale_axial)
            if shear:
                el.plot_shear_force(ax=ax, analysis_case=analysis_case, scalef=scale_shear)
            if moment:
                el.plot_bending_moment(ax=ax, analysis_case=analysis_case, scalef=scale_moment)

        # plot the undeformed structure
        self.plot_geom(analysis_case=analysis_case, ax=ax)

    def plot_buckling_results(self, analysis_case, buckling_mode=0):
        """Method used to plot a buckling eigenvector. The undeformed structure is plotted with a
        dashed line.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int buckling_mode: Buckling mode to plot
        """

        (fig, ax) = plt.subplots()

        # set initial plot limits
        (xmin, xmax, ymin, ymax, _, _) = self.analysis.get_node_lims()

        # determine max eigenvector displacement value (ignore rotation)
        max_v = 0

        # loop through all the elements
        for el in self.analysis.elements:
            (w, v_el) = el.get_buckling_results(
                analysis_case=analysis_case, buckling_mode=buckling_mode)

            max_v = max(max_v, abs(v_el[0, 0]), abs(v_el[0, 1]), abs(v_el[1, 0]), abs(v_el[1, 1]))

        # determine plot scale
        scale = 0.1 * max(xmax - xmin, ymax - ymin) / max_v

        # plot eigenvectors
        for el in self.analysis.elements:
            (_, v_el) = el.get_buckling_results(
                analysis_case=analysis_case, buckling_mode=buckling_mode)

            el.plot_deformed_element(
                ax=ax, analysis_case=analysis_case, n=self.n_subdiv, def_scale=scale, u_el=v_el)

        # plot the load factor (eigenvalue)
        ax.set_title("Load Factor for Mode {:d}: {:.4e}".format(buckling_mode, w), size=10)

        # plot the undeformed structure
        self.plot_geom(analysis_case=analysis_case, ax=ax, dashed=True)

    def plot_frequency_results(self, analysis_case, frequency_mode=0):
        """Method used to plot a frequency eigenvector. The undeformed structure is plotted with a
        dashed line.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param int frequency_mode: Frequency mode to plot
        """

        (fig, ax) = plt.subplots()

        # set initial plot limits
        (xmin, xmax, ymin, ymax, _, _) = self.analysis.get_node_lims()

        # determine max eigenvector displacement value (ignore rotation)
        max_v = 0

        # loop through all the elements
        for el in self.analysis.elements:
            (w, v_el) = el.get_frequency_results(
                analysis_case=analysis_case, frequency_mode=frequency_mode)

            max_v = max(max_v, abs(v_el[0, 0]), abs(v_el[0, 1]), abs(v_el[1, 0]), abs(v_el[1, 1]))

        # determine plot scale
        scale = 0.1 * max(xmax - xmin, ymax - ymin) / max_v

        # plot eigenvectors
        for el in self.analysis.elements:
            (_, v_el) = el.get_frequency_results(
                analysis_case=analysis_case, frequency_mode=frequency_mode)

            el.plot_deformed_element(
                ax=ax, analysis_case=analysis_case, n=self.n_subdiv, def_scale=scale, u_el=v_el)

        # plot the frequency (eigenvalue)
        ax.set_title("Frequency for Mode {:d}: {:.4e} Hz".format(frequency_mode, w), size=10)

        # plot the undeformed structure
        self.plot_geom(analysis_case=analysis_case, ax=ax, dashed=True)

    def get_support_angle(self, node, prefer_dir=None):
        """Given a node object, returns the optimal angle to plot a support. Essentially finds the
        average angle of the connected elements and considers a preferred plotting direction. N.B.
        this method is adopted from the MATLAB code by F.P. van der Meer: plotGeom.m.

        :param node: Node object
        :type node: :class:`~feastruct.fea.node.node`
        :param int prefer_dir: Preferred direction to plot the support, where 0 corresponds to the
            x-axis and 1 corresponds to the y-axis
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
            if prefer_dir == 1:
                if max(np.sin([phi[i0] * np.pi / 180, phi[i0+1] * np.pi / 180])) > -0.1:
                    angle = 90
            elif prefer_dir == 0:
                if max(np.cos([phi[i0] * np.pi / 180, phi[i0+1] * np.pi / 180])) > -0.1:
                    angle = 0

        return (angle, num_el)

    def wide_lim(self, x):
        """Returns a tuple corresponding to the axis limits (x) stretched by 2% on either side.

        :param x: List containing axis limits e.g. [xmin, xmax]
        :type x: list[float, float]
        :return: Stretched axis limits (x1, x2)
        :rtype: tuple(float, float)
        """

        x2 = max(x)
        x1 = min(x)
        dx = x2-x1
        f = 0.02

        return (x1 - f * dx, x2 + f * dx)


class ScalarResult:
    """Class for storing a scalar result for a specific analysis case.

    :cvar float result: Scalar result
    :cvar analysis_case: Analysis case relating to the result
    :vartype analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
    """

    def __init__(self, result, analysis_case):
        """Inits the ScalarResult class.

        :cvar float result: Scalar result
        :param analysis_case: Analysis case relating to the result
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        self.result = result
        self.analysis_case = analysis_case
