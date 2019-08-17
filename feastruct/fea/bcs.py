import numpy as np
from matplotlib.patches import Polygon
from feastruct.post.post import ScalarResult


class BoundaryCondition:
    """Parent class for supports and loads.

    Provides an init method for the creation of boundary conditions.

    :cvar node: The node object at which the boundary condition acts
    :vartype node: :class:`~feastruct.fea.node.Node`
    :cvar float val: The value of the boundary condition
    :cvar int dof: The degree of freedom about which the boundary condition acts
    """

    def __init__(self, node, val, dof):
        """Inits the BoundaryCondition class.

        :param node: The node object at which the boundary condition acts
        :type node: :class:`~feastruct.fea.node.Node`
        :param float val: The value of the boundary condition
        :param int dof: The degree of freedom about which the boundary condition acts
        """

        # assign the node object, value and dof of the boundary condition
        self.node = node
        self.val = val
        self.dof = dof

    def get_gdof(self):
        """Returns the global degree of freedom number for the boundary condition.

        :returns: Global degree of freedom number
        :rtype: int
        """

        return self.node.dofs[self.dof].global_dof_num


class NodalSupport(BoundaryCondition):
    """Class for a dirichlet boundary condition acting at a node.

    Provides methods for the FEA solver and post-processing.

    :cvar node: The node object at which the boundary condition acts
    :vartype node: :class:`~feastruct.fea.node.Node`
    :cvar float val: The value of the boundary condition
    :cvar int dof: The degree of freedom about which the boundary condition acts
    :cvar reactions: A list of reaction objects
    :vartype reactions: list[:class:`~feastruct.post.post.ScalarResult`]
    """

    def __init__(self, node, val, dof):
        """inits the NodalSupport class.

        :param node: The node object at which the boundary condition acts
        :type node: :class:`~feastruct.fea.node.Node`
        :param float val: The value of the boundary condition
        :param int dof: The degree of freedom about which the boundary condition acts
        """

        # initialise the parent class
        super().__init__(node, val, dof)

        # initialise the nodal reaction results
        self.reactions = []

    def apply_support(self, K, f_ext):
        """Applies the nodal support.

        The stiffness matrix and external force vector are modified to apply the dirichlet boundary
        condition to enforce the displacement at the chosen degree of freedom to be equal to the
        specified value.

        :param K: Global stiffness matrix of size *(N x N)*
        :type K: :class:`numpy.ndarray`
        :param f_ext: Global external force vector of size *N*
        :type f_ext: :class:`numpy.ndarray`
        """

        # get gdof number for the support
        gdof = self.node.dofs[self.dof].global_dof_num

        # modify stiffness matrix and f_ext
        K[gdof, :] = 0
        K[gdof, gdof] = 1
        f_ext[gdof] = self.val

    def get_reaction(self, analysis_case):
        """Gets the reaction force result corresponding to analysis_case.

        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :returns: Reaction force at the node
        :rtype: float
        """

        # loop through reactions
        for reaction in self.reactions:
            if reaction.analysis_case == analysis_case:
                return reaction.result

    def save_reaction(self, f, analysis_case):
        """Saves the reaction force corresponding to analysis_case.

        :param float f: Reaction force at the node
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # check to see if there is already a reaction for the current analysis_case
        for reaction in self.reactions:
            if reaction.analysis_case == analysis_case:
                reaction.result = f
                return

        # if there isn't already a reaction for the current analysis_case
        self.reactions.append(ScalarResult(result=f, analysis_case=analysis_case))

    def plot_support(self, ax, small, get_support_angle, analysis_case, deformed, def_scale):
        """Plots a graphical representation of the nodal support.

        Based on the type of support at the node, a graphical representation of the support type is
        generated and plotted. Possible support types include rollers, hinges, rotation restraints,
        fixed rollers and fully fixed supports. The angle of the connecting elements is considered
        in order to produce the most visually appealing representation. The support location is
        displaced if a deformed plot is desired. Note that some of the methods used to plot the
        supports are taken from Frans van der Meer's code plotGeom.m.

        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param float small: A dimension used to scale the support
        :param get_support_angle: A function that returns the support angle and the number of
            connected elements
        :type get_support_angle: :func:`feastruct.post.post.PostProcessor.get_support_angle`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param bool deformed: Represents whether or not the node locations are deformed based on
            the results of analysis_case
        :param float def_scale: Value used to scale deformations
        """

        fixity = analysis_case.freedom_case.get_nodal_fixities(node=self.node)

        if fixity not in ([1, 1, 0], [0, 1, 0]):
            (angle, num_el) = get_support_angle(self.node)

        if fixity == [1, 0, 0]:
            # ploy a y-roller
            angle = round(angle / 180) * 180
            self.plot_xysupport(ax, angle, True, num_el == 1, small, analysis_case, deformed,
                                def_scale)

        elif fixity == [0, 1, 0]:
            (angle, num_el) = get_support_angle(self.node, 1)
            # plot an x-roller
            if np.mod(angle + 1, 180) < 2:  # prefer support below
                angle = 90
            else:
                angle = round((angle + 90) / 180) * 180 - 90

            self.plot_xysupport(ax, angle, True, num_el == 1, small, analysis_case, deformed,
                                def_scale)

        elif fixity == [1, 1, 0]:
            # plot a hinge
            (angle, num_el) = get_support_angle(self.node, 1)
            self.plot_xysupport(ax, angle, False, num_el == 1, small, analysis_case, deformed,
                                def_scale)

        elif fixity == [0, 0, 1]:
            ax.plot(self.node.x, self.node.y, 'kx', markersize=8)

        else:
            # plot a support with moment fixity
            if fixity == [1, 1, 1]:
                # plot a fixed support
                s = np.sin(angle * np.pi / 180)
                c = np.cos(angle * np.pi / 180)
                rot_mat = np.array([[c, -s], [s, c]])
                line = np.array([[0, 0], [-1, 1]]) * small
                rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * small
                ec = 'none'

            elif fixity == [1, 0, 1]:
                # plot y-roller block
                angle = round(angle / 180) * 180
                s = np.sin(angle * np.pi / 180)
                c = np.cos(angle * np.pi / 180)
                rot_mat = np.array([[c, -s], [s, c]])
                line = np.array([[-0.85, -0.85], [-1, 1]]) * small
                rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * small
                ec = 'k'

            elif fixity == [0, 1, 1]:
                # plot x-roller block
                angle = round((angle + 90) / 180) * 180 - 90
                s = np.sin(angle * np.pi / 180)
                c = np.cos(angle * np.pi / 180)
                rot_mat = np.array([[c, -s], [s, c]])
                line = np.array([[-0.85, -0.85], [-1, 1]]) * small
                rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * small
                ec = 'k'

            rot_line = np.matmul(rot_mat, line)
            rot_rect = np.matmul(rot_mat, rect)

            # add coordinates of node
            if deformed:
                # get displacement of node for current analysis case
                u = [0, 0]
                u[0] = self.node.dofs[0].get_displacement(analysis_case)
                u[1] = self.node.dofs[1].get_displacement(analysis_case)

                rot_line[0, :] += self.node.x + u[0] * def_scale
                rot_line[1, :] += self.node.y + u[1] * def_scale
                rot_rect[0, :] += self.node.x + u[0] * def_scale
                rot_rect[1, :] += self.node.y + u[1] * def_scale
            else:
                rot_line[0, :] += self.node.x
                rot_line[1, :] += self.node.y
                rot_rect[0, :] += self.node.x
                rot_rect[1, :] += self.node.y

            ax.plot(rot_line[0, :], rot_line[1, :], 'k-', linewidth=1)
            ax.add_patch(Polygon(np.transpose(rot_rect), facecolor=(0.7, 0.7, 0.7), edgecolor=ec))

    def plot_imposed_disp(self, ax, max_disp, small, get_support_angle, analysis_case, deformed,
                          def_scale):
        """Plots a graphical representation of an imposed translation.

        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param float max_disp: Maximum imposed displacement in the analysis case
        :param float small: A dimension used to scale the support
        :param get_support_angle: A function that returns the support angle and the number of
            connected elements
        :type get_support_angle: :func:`feastruct.post.post.PostProcessor.get_support_angle`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param bool deformed: Represents whether or not the node locations are deformed based on
            the results of case id
        :param float def_scale: Value used to scale deformations
        """

        val = self.val / max_disp
        offset = 0.5 * small

        lf = abs(val) * 1.5 * small  # arrow length
        lh = 0.6 * small  # arrow head length
        wh = 0.6 * small  # arrow head width
        sp = 0.15 * small  # half spacing between double line
        lf = max(lf, lh * 1.5)

        (angle, num_el) = get_support_angle(self.node)
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        n = np.array([c, s])
        inward = (n[self.dof] == 0 or np.sign(n[self.dof]) == np.sign(val))

        to_rotate = (self.dof) * 90 + (n[self.dof] >= 0) * 180
        sr = np.sin(to_rotate * np.pi / 180)
        cr = np.cos(to_rotate * np.pi / 180)
        rot_mat = np.array([[cr, -sr], [sr, cr]])

        x0 = offset + inward * lf
        x2 = offset + (not inward) * lf
        x1 = x2 + (inward) * lh - (not inward) * lh
        pp = np.array([[x1, x1, x2], [-wh / 2, wh / 2, 0]])
        ll = np.array([[x1, x0, x0, x1], [sp, sp, -sp, -sp]])

        rl = np.matmul(rot_mat, ll)
        rp = np.matmul(rot_mat, pp)

        # add coordinates of node
        if deformed:
            # get displacement of node for current analysis case
            u = [0, 0]
            u[0] = self.node.dofs[0].get_displacement(analysis_case)
            u[1] = self.node.dofs[1].get_displacement(analysis_case)

            rp[0, :] += self.node.x + u[0] * def_scale
            rp[1, :] += self.node.y + u[1] * def_scale
            rl[0, :] += self.node.x + u[0] * def_scale
            rl[1, :] += self.node.y + u[1] * def_scale
        else:
            rp[0, :] += self.node.x
            rp[1, :] += self.node.y
            rl[0, :] += self.node.x
            rl[1, :] += self.node.y

        ax.plot(rl[0, :], rl[1, :], 'k-')
        ax.add_patch(Polygon(np.transpose(rp), facecolor='none', linewidth=1, edgecolor='k'))

    def plot_imposed_rot(self, ax, small, get_support_angle, analysis_case, deformed, def_scale):
        """Plots a graphical representation of an imposed rotation.

        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param float small: A dimension used to scale the support
        :param get_support_angle: A function that returns the support angle and the number of
            connected elements
        :type get_support_angle: :func:`feastruct.post.post.PostProcessor.get_support_angle`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param bool deformed: Represents whether or not the node locations are deformed based on
            the results of case id
        :param float def_scale: Value used to scale deformations
        """

        lh = 0.4 * small  # arrow head length
        wh = 0.4 * small  # arrow head width
        r1 = 1.0 * small
        r2 = 1.2 * small
        (angle, num_el) = get_support_angle(self.node)
        ths = np.arange(100, 261)

        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])

        # make arrow tail around (0,0)
        rr = (r1 + r2) / 2
        ll = np.array([rr * np.cos(ths * np.pi / 180), rr * np.sin(ths * np.pi / 180)])
        l1 = np.array([r1 * np.cos(ths * np.pi / 180), r1 * np.sin(ths * np.pi / 180)])
        l2 = np.array([r2 * np.cos(ths * np.pi / 180), r2 * np.sin(ths * np.pi / 180)])

        # make arrow head at (0,0)
        pp = np.array([[-lh, -lh, 0], [-wh / 2, wh / 2, 0]])

        # rotate arrow head around (0,0)
        if self.val > 0:
            thTip = 90 - ths[11]
            xTip = ll[:, -1]
            l1 = l1[:, 1:-21]
            l2 = l2[:, 1:-21]
            ibase = 0
        else:
            thTip = ths[11] - 90
            xTip = ll[:, 0]
            l1 = l1[:, 21:]
            l2 = l2[:, 21:]
            ibase = np.shape(l1)[1] - 1

        cTip = np.cos(thTip * np.pi / 180)
        sTip = np.sin(thTip * np.pi / 180)
        rTip = np.array([[cTip, -sTip], [sTip, cTip]])
        pp = np.matmul(rTip, pp)

        # shift arrow head to tip
        pp[0, :] += xTip[0]
        pp[1, :] += xTip[1]

        # rotate arrow to align it with the node
        rl1 = np.matmul(rot_mat, l1)
        rl2 = np.matmul(rot_mat, l2)
        rp = np.matmul(rot_mat, pp)

        # add coordinates of node
        if deformed:
            # get displacement of node for current analysis case
            u = [0, 0]
            u[0] = self.node.dofs[0].get_displacement(analysis_case)
            u[1] = self.node.dofs[1].get_displacement(analysis_case)

            rp[0, :] += self.node.x + u[0] * def_scale
            rp[1, :] += self.node.y + u[1] * def_scale
            rl1[0, :] += self.node.x + u[0] * def_scale
            rl1[1, :] += self.node.y + u[1] * def_scale
            rl2[0, :] += self.node.x + u[0] * def_scale
            rl2[1, :] += self.node.y + u[1] * def_scale
        else:
            rp[0, :] += self.node.x
            rp[1, :] += self.node.y
            rl1[0, :] += self.node.x
            rl1[1, :] += self.node.y
            rl2[0, :] += self.node.x
            rl2[1, :] += self.node.y

        # shift arrow to node and plot
        ax.plot(rl1[0, :], rl1[1, :], 'k-')
        ax.plot(rl2[0, :], rl2[1, :], 'k-')
        ax.plot(np.append(rl1[0, ibase], rl2[0, ibase]), np.append(rl1[1, ibase], rl2[1, ibase]),
                'k-')
        ax.add_patch(Polygon(np.transpose(rp), facecolor='none', linewidth=1, edgecolor='k'))

    def plot_reaction(self, ax, max_reaction, small, get_support_angle, analysis_case):
        """Plots a graphical representation of a reaction force and displays the value of the
        reaction force. A straight arrow is plotted for a translational reaction and a curved arrow
        is plotted for a rotational reaction.

        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param float max_reaction: Maximum reaction force in the analysis case
        :param float small: A dimension used to scale the support
        :param get_support_angle: A function that returns the support angle and the number of
            connected elements
        :type get_support_angle: :func:`feastruct.post.post.PostProcessor.get_support_angle`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        """

        # get reaction force
        reaction = self.get_reaction(analysis_case)

        # dont plot small reaction
        if abs(reaction) < 1e-6:
            return

        if self.dof in [0, 1]:
            val = reaction / max_reaction

            lf = abs(val) * 1.5 * small  # arrow length
            lh = 0.4 * small  # arrow head length
            wh = 0.4 * small  # arrow head width
            lf = max(lf, lh * 1.5)
            offset = 0.5 * small
            xoff = 0
            yoff = 0

            if self.dof == 0:
                rot_mat = np.array([[-1, 0], [0, -1]]) * np.sign(val)
                va = 'center'
                if val > 0:
                    ha = 'right'
                    xoff = -offset / 2
                else:
                    ha = 'left'
                    xoff = offset / 2
            elif self.dof == 1:
                rot_mat = np.array([[0, 1], [-1, 0]]) * np.sign(val)
                ha = 'center'
                if val > 0:
                    va = 'top'
                    yoff = -offset / 2
                else:
                    va = 'bottom'
                    yoff = offset / 2

            inward = True

            ll = np.array([[offset, offset + lf], [0, 0]])
            p0 = offset
            p1 = p0 + lh
            pp = np.array([[p1, p1, p0], [-wh / 2, wh / 2, 0]])

            # correct end of arrow line
            if inward:
                ll[0, 0] += lh
            else:
                ll[0, 1] -= lh

            rl = np.matmul(rot_mat, ll)
            rp = np.matmul(rot_mat, pp)
            rl[0, :] += self.node.x
            rl[1, :] += self.node.y
            rp[0, :] += self.node.x
            rp[1, :] += self.node.y
            s = 0
            e = None

            tl = np.array([rl[0, 1] + xoff, rl[1, 1] + yoff])
        else:
            (angle, num_el) = get_support_angle(self.node)
            s = np.sin(angle * np.pi / 180)
            c = np.cos(angle * np.pi / 180)

            lh = 0.3 * small  # arrow head length
            wh = 0.3 * small  # arrow head width
            rr = 1.5 * small
            ths = np.arange(100, 261)
            rot_mat = np.array([[c, -s], [s, c]])

            # make arrow tail around (0,0)
            ll = np.array([rr * np.cos(ths * np.pi / 180),
                           rr * np.sin(ths * np.pi / 180)])

            # make arrow head at (0,0)
            pp = np.array([[-lh, -lh, 0], [-wh / 2, wh / 2, 0]])

            # rotate arrow head around (0,0)
            if reaction > 0:
                thTip = 90 - ths[11]
                xTip = ll[:, -1]
                s = 0
                e = -1
            else:
                thTip = ths[11] - 90
                xTip = ll[:, 0]
                s = 1
                e = None

            cTip = np.cos(thTip * np.pi / 180)
            sTip = np.sin(thTip * np.pi / 180)
            rTip = np.array([[cTip, -sTip], [sTip, cTip]])
            pp = np.matmul(rTip, pp)

            # shift arrow head to tip
            pp[0, :] += xTip[0]
            pp[1, :] += xTip[1]

            # rotate arrow to align it with the node
            rl = np.matmul(rot_mat, ll)
            rp = np.matmul(rot_mat, pp)
            rl[0, :] += self.node.x
            rl[1, :] += self.node.y
            rp[0, :] += self.node.x
            rp[1, :] += self.node.y

            ha = 'center'
            va = 'center'
            tl = np.array([rl[0, 0], rl[1, 0]])

        ax.plot(rl[0, s:e], rl[1, s:e], linewidth=1.5, color='r')
        ax.add_patch(Polygon(np.transpose(rp), facecolor='r'))
        ax.text(tl[0], tl[1], "{:5.3g}".format(reaction), size=8, horizontalalignment=ha,
                verticalalignment=va)

    def plot_xysupport(self, ax, angle, roller, hinge, small, analysis_case, deformed, def_scale):
        """Plots a hinged or roller support.

        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param float angle: Angle at which to plot the support
        :param bool roller: Whether or not the support is a roller
        :param bool hinge: Whether or not there is a hinge at the support
        :param float small: A dimension used to scale the support
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param bool deformed: Represents whether or not the node locations are deformed based on
            the results of case id
        :param float def_scale: Value used to scale deformations
        """

        # determine coordinates of node
        if deformed:
            # get displacement of node for current analysis case
            u = [0, 0]
            u[0] = self.node.dofs[0].get_displacement(analysis_case)
            u[1] = self.node.dofs[1].get_displacement(analysis_case)

            x = self.node.x + u[0] * def_scale
            y = self.node.y + u[1] * def_scale
        else:
            x = self.node.x
            y = self.node.y

        # determine coordinates of triangle
        dx = small
        h = np.sqrt(3) / 2
        triangle = np.array([[-h, -h, -h, 0, -h], [-1, 1, 0.5, 0, -0.5]]) * dx
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])
        rot_triangle = np.matmul(rot_mat, triangle)

        if roller:
            line = np.array([[-1.1, -1.1], [-1, 1]]) * dx
            rot_line = np.matmul(rot_mat, line)
            ax.plot(rot_line[0, :] + x, rot_line[1, :] + y, 'k-', linewidth=1)
        else:
            rect = np.array([[-1.4, -1.4, -h, -h], [-1, 1, 1, -1]]) * dx
            rot_rect = np.matmul(rot_mat, rect)
            rot_rect[0, :] += x
            rot_rect[1, :] += y
            ax.add_patch(Polygon(np.transpose(rot_rect), facecolor=(0.7, 0.7, 0.7)))

        ax.plot(rot_triangle[0, :] + x, rot_triangle[1, :] + y, 'k-', linewidth=1)

        if hinge:
            ax.plot(x, y, 'ko', markerfacecolor='w', linewidth=1, markersize=4)


class NodalLoad(BoundaryCondition):
    """Class for a neumann boundary condition acting at a node.

    Provides methods for the FEA solver and post-processing.

    :cvar node: The node object at which the nodal load acts
    :vartype node: :class:`~feastruct.fea.node.Node`
    :cvar float val: The value of the nodal load
    :cvar int dof: The degree of freedom about which the nodal load acts
    """

    def __init__(self, node, val, dof):
        """inits the NodalLoad class.

        :param node: The node object at which the nodal load acts
        :type node: :class:`~feastruct.fea.node.Node`
        :param float val: The value of the nodal load
        :param int dof: The degree of freedom about which the nodal load acts
        """

        super().__init__(node, val, dof)

    def apply_load(self, f_ext):
        """Applies the nodal load.

        The external force vector is modified to apply the neumann boundary condition.

        :param f_ext: Global external force vector of size *N*
        :type f_ext: :class:`numpy.ndarray`
        """

        # get gdof number for the load
        gdof = self.node.dofs[self.dof].global_dof_num

        # add load to f_ext, selecting the correct dof from dofs
        f_ext[gdof] += self.val

    def plot_load(self, ax, max_force, small, get_support_angle, analysis_case, deformed,
                  def_scale):
        """Plots a graphical representation of a nodal force. A straight arrow is plotted for a
        translational load and a curved arrow is plotted for a moment.

        :param ax: Axes object on which to plot
        :type ax: :class:`matplotlib.axes.Axes`
        :param float max_force: Maximum translational nodal load
        :param float small: A dimension used to scale the support
        :param get_support_angle: A function that returns the support angle and the number of
            connected elements
        :type get_support_angle: :func:`feastruct.post.post.PostProcessor.get_support_angle`
        :param analysis_case: Analysis case
        :type analysis_case: :class:`~feastruct.fea.cases.AnalysisCase`
        :param bool deformed: Represents whether or not the node locations are deformed based on
            the results of case id
        :param float def_scale: Value used to scale deformations
        """

        # determine coordinates of node
        if deformed:
            u = [0, 0]
            u[0] = self.node.dofs[0].get_displacement(analysis_case)
            u[1] = self.node.dofs[1].get_displacement(analysis_case)

            x = self.node.x + u[0] * def_scale
            y = self.node.y + u[1] * def_scale
        else:
            x = self.node.x
            y = self.node.y

        if max_force == 0:
            val = self.val
        else:
            val = self.val / max_force

        offset = 0.5 * small
        (angle, num_el) = get_support_angle(self.node)
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)

        # plot nodal force
        if self.dof in [0, 1]:
            lf = abs(val) * 1.5 * small  # arrow length
            lh = 0.6 * small  # arrow head length
            wh = 0.6 * small  # arrow head width
            lf = max(lf, lh * 1.5)

            n = np.array([c, s])
            inward = (n[self.dof] == 0 or np.sign(n[self.dof]) == np.sign(val))

            to_rotate = (self.dof) * 90 + (n[self.dof] > 0) * 180
            sr = np.sin(to_rotate * np.pi / 180)
            cr = np.cos(to_rotate * np.pi / 180)
            rot_mat = np.array([[cr, -sr], [sr, cr]])

            ll = np.array([[offset, offset + lf], [0, 0]])
            p0 = offset + (not inward) * lf
            p1 = p0 + (inward) * lh - (not inward) * lh
            pp = np.array([[p1, p1, p0], [-wh / 2, wh / 2, 0]])

            # correct end of arrow line
            if inward:
                ll[0, 0] += lh
            else:
                ll[0, 1] -= lh

            rl = np.matmul(rot_mat, ll)
            rp = np.matmul(rot_mat, pp)
            rp[0, :] += x
            rp[1, :] += y
            s = 0
            e = None

        # plot nodal moment
        else:
            lh = 0.4 * small  # arrow head length
            wh = 0.4 * small  # arrow head width
            rr = 1.5 * small
            ths = np.arange(100, 261)
            rot_mat = np.array([[c, -s], [s, c]])

            # make arrow tail around (0,0)
            ll = np.array([rr * np.cos(ths * np.pi / 180), rr * np.sin(ths * np.pi / 180)])

            # make arrow head at (0,0)
            pp = np.array([[-lh, -lh, 0], [-wh / 2, wh / 2, 0]])

            # rotate arrow head around (0,0)
            if val > 0:
                thTip = 90 - ths[11]
                xTip = ll[:, -1]
                s = 0
                e = -1
            else:
                thTip = ths[11] - 90
                xTip = ll[:, 0]
                s = 1
                e = None

            cTip = np.cos(thTip * np.pi / 180)
            sTip = np.sin(thTip * np.pi / 180)
            rTip = np.array([[cTip, -sTip], [sTip, cTip]])
            pp = np.matmul(rTip, pp)

            # shift arrow head to tip
            pp[0, :] += xTip[0]
            pp[1, :] += xTip[1]

            # rotate arrow to align it with the node
            rl = np.matmul(rot_mat, ll)
            rp = np.matmul(rot_mat, pp)
            rp[0, :] += x
            rp[1, :] += y

        ax.plot(rl[0, s:e] + x, rl[1, s:e] + y, 'k-', linewidth=2)
        ax.add_patch(Polygon(np.transpose(rp), facecolor='k'))


class ElementLoad:
    """Parent class for a load applied to a finite element.

    :cvar element: Finite element to which the load is applied
    :vartype element: :class:`~feastruct.fea.fea.FiniteElement`
    """

    def __init__(self, element):
        """Inits the ElementLoad class.

        :param element: Finite element to which the load is applied
        :type element: :class:`~feastruct.fea.fea.FiniteElement`
        """

        self.element = element

    def nodal_equivalent_loads(self):
        """Placehodler for the nodal_equivalent_loads method."""

        pass

    def apply_load(self, f_eq):
        """Placeholder for the apply_load method.

        :param f_eq: Global equivalent nodal loads vector of size *N*
        :type f_eq: :class:`numpy.ndarray
        """

        pass

    def plot_load(self):
        """Placeholder for the plot_load method."""

        pass
