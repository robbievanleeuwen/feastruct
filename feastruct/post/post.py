import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import pdb


class PostProcessor:
    """docstring for PostProcessor.
    """

    def __init__(self, analysis):
        """asldkjasld
        """

        self.analysis = analysis

    def plot_geom(self):
        """askldjasld

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        (fig, ax) = plt.subplots()

        # TODO: plot deformation

        for el in self.analysis.elements:
            el.plot_element(ax)

        # set initial plot limits
        (xmin, xmax, ymin, ymax) = self.analysis.get_node_lims()
        ax.set_xlim(xmin-1e-12, xmax)
        ax.set_ylim(ymin-1e-12, ymax)

        # get 2% of the maxmimum dimension
        self.get_small = 0.02*max(xmax-xmin, ymax-ymin)

        # generate lists of nodal supports and imposed displacements
        support_node_list = []
        imposed_disp_list = []
        max_disp = 0

        for support in self.analysis.supports:
            # check that there is no imposed displacement
            if support["val"] == 0:
                support["node"].fixity[support["dir"]-1] = 1

                if support["node"] not in support_node_list:
                    support_node_list.append(support["node"])
            # if there is an imposed displacement
            else:
                imposed_disp_list.append(support)
                if support["dir"] == 1 or support["dir"] == 2:
                    max_disp = max(max_disp, abs(support["val"]))

        # plot supports
        for node in support_node_list:
            if node.fixity == [1, 0, 0]:
                self.plot_yroller(ax, node)
            if node.fixity == [0, 1, 0]:
                self.plot_xroller(ax, node)
            if node.fixity == [1, 1, 0]:
                self.plot_hinge(ax, node)
            if node.fixity == [1, 1, 1]:
                self.plot_fixed(ax, node)
            if node.fixity == [1, 0, 1]:
                self.plot_yroller_block(ax, node)
            if node.fixity == [0, 1, 1]:
                self.plot_xroller_block(ax, node)
            if node.fixity == [0, 0, 1]:
                self.plot_zmoment(ax, node)

        # plot imposed displacements
        for imposed_disp in imposed_disp_list:
            if imposed_disp["dir"] == 1 or imposed_disp["dir"] == 2:
                self.plot_imposed_disp(ax, imposed_disp, max_disp)
            elif imposed_disp["dir"] == 3:
                self.plot_imposed_rot(ax, imposed_disp)

        # find max nodal loads
        max_force = 0

        for nodal_load in self.analysis.nodal_loads:
            if nodal_load["dir"] == 1 or nodal_load["dir"] == 2:
                max_force = max(max_force, abs(nodal_load["val"]))

        # plot nodal loads
        for nodal_load in self.analysis.nodal_loads:
            if nodal_load["dir"] == 1 or nodal_load["dir"] == 2:
                self.plot_nodal_force(ax, nodal_load, max_force)
            elif nodal_load["dir"] == 3:
                self.plot_nodal_moment(ax, nodal_load)

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

    def plot_xroller(self, ax, node):
        """asdksa

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        (angle, num_el) = self.get_support_angle(node)

        # prefer support below
        if np.mod(angle + 1, 180) < 2:
            angle = 90
        else:
            angle = round((angle + 90) / 180) * 180 - 90

        self.plot_xysupport(ax, node, angle, True, num_el == 1)

    def plot_yroller(self, ax, node):
        """asdksa

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        (angle, num_el) = self.get_support_angle(node)
        angle = round(angle / 180) * 180
        self.plot_xysupport(ax, node, angle, True, num_el == 1)

    def plot_hinge(self, ax, node):
        """asdksa

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        (angle, num_el) = self.get_support_angle(node, 2)
        self.plot_xysupport(ax, node, angle, False, num_el == 1)

    def plot_fixed(self, ax, node):
        """slkdjaksldasd
        """

        dx = self.get_small
        (angle, num_el) = self.get_support_angle(node)
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])
        line = np.array([[0, 0], [-1, 1]]) * dx
        rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * dx
        rot_line = np.matmul(rot_mat, line)
        rot_rect = np.matmul(rot_mat, rect)
        rot_rect[0, :] += node.x
        rot_rect[1, :] += node.y

        ax.plot(rot_line[0, :] + node.x, rot_line[1, :] + node.y, 'k-',
                linewidth=1)
        ax.add_patch(Polygon(np.transpose(rot_rect),
                             facecolor=(0.7, 0.7, 0.7)))

    def plot_xroller_block(self, ax, node):
        """slkdjaksldasd
        """

        dx = self.get_small
        (angle, num_el) = self.get_support_angle(node)
        angle = round((angle + 90) / 180) * 180 - 90
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])
        line = np.array([[-0.85, -0.85], [-1, 1]]) * dx
        rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * dx
        rot_line = np.matmul(rot_mat, line)
        rot_rect = np.matmul(rot_mat, rect)
        rot_rect[0, :] += node.x
        rot_rect[1, :] += node.y

        ax.plot(rot_line[0, :] + node.x, rot_line[1, :] + node.y, 'k-',
                linewidth=1)
        ax.add_patch(Polygon(np.transpose(rot_rect),
                             facecolor=(0.7, 0.7, 0.7), edgecolor='k'))

    def plot_yroller_block(self, ax, node):
        """slkdjaksldasd
        """

        dx = self.get_small
        (angle, num_el) = self.get_support_angle(node)
        angle = round(angle / 180) * 180
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])
        line = np.array([[-0.85, -0.85], [-1, 1]]) * dx
        rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * dx
        rot_line = np.matmul(rot_mat, line)
        rot_rect = np.matmul(rot_mat, rect)
        rot_rect[0, :] += node.x
        rot_rect[1, :] += node.y

        ax.plot(rot_line[0, :] + node.x, rot_line[1, :] + node.y, 'k-',
                linewidth=1)
        ax.add_patch(Polygon(np.transpose(rot_rect),
                             facecolor=(0.7, 0.7, 0.7), edgecolor='k'))

    def plot_zmoment(self, ax, node):
        """slkdjaksldasd
        """

        ax.plot(node.x, node.y, 'kx', markersize=8)

    def plot_xysupport(self, ax, node, angle, roller, hinge):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        # determine coordinates of triangle
        dx = self.get_small
        h = np.sqrt(3) / 2
        triangle = np.array([[-h, -h, -h, 0, -h], [-1, 1, 0.5, 0, -0.5]]) * dx
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])
        rot_triangle = np.matmul(rot_mat, triangle)

        if roller:
            line = np.array([[-1.1, -1.1], [-1, 1]]) * dx
            rot_line = np.matmul(rot_mat, line)
            ax.plot(rot_line[0, :] + node.x, rot_line[1, :] + node.y, 'k-',
                    linewidth=1)
        else:
            rect = np.array([[-1.4, -1.4, -h, -h], [-1, 1, 1, -1]]) * dx
            rot_rect = np.matmul(rot_mat, rect)
            rot_rect[0, :] += node.x
            rot_rect[1, :] += node.y
            ax.add_patch(Polygon(np.transpose(rot_rect),
                                 facecolor=(0.7, 0.7, 0.7)))

        ax.plot(rot_triangle[0, :] + node.x, rot_triangle[1, :] + node.y, 'k-',
                linewidth=1)

        if hinge:
            ax.plot(node.x, node.y, 'ko', markerfacecolor='w',
                    linewidth=1, markersize=4)

    def plot_nodal_force(self, ax, nodal_load, max_force):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        val = nodal_load["val"] / max_force
        node = nodal_load["node"]
        dir = nodal_load["dir"]

        small = self.get_small
        offset = 0.5 * small

        lf = abs(val) * 1.5 * small  # arrow length
        lh = 0.6 * small  # arrow head length
        wh = 0.6 * small  # arrow head width
        lf = max(lf, lh * 1.5)

        (angle, num_el) = self.get_support_angle(node)
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        n = np.array([c, s])
        inward = (n[dir-1] == 0 or np.sign(n[dir-1]) == np.sign(val))

        to_rotate = (dir - 1) * 90 + (n[dir-1] > 0) * 180
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
        rp[0, :] += node.x
        rp[1, :] += node.y

        ax.plot(rl[0, :] + node.x, rl[1, :] + node.y, 'k-', linewidth=2)
        ax.add_patch(Polygon(np.transpose(rp), facecolor='k'))

    def plot_nodal_moment(self, ax, nodal_load):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        node = nodal_load["node"]
        val = nodal_load["val"]

        small = self.get_small
        lh = 0.4 * small  # arrow head length
        wh = 0.4 * small  # arrow head width
        rr = 1.5 * small
        (angle, num_el) = self.get_support_angle(node)
        ths = np.arange(100, 261)

        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])

        # make arrow tail around (0,0)
        ll = np.array([rr * np.cos(ths * np.pi / 180),
                       rr * np.sin(ths * np.pi / 180)])

        # make arrow head at (0,0)
        pp = np.array([[-lh, -lh, 0], [-wh / 2, wh / 2, 0]])

        # rotate arrow head around (0,0)
        if val > 0:
            thTip = 90 - ths[11]
            xTip = ll[:, -1]
            s = 0
            e = -2
        else:
            thTip = ths[11] - 90
            xTip = ll[:, 0]
            s = 1
            e = -1

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
        rp[0, :] += node.x
        rp[1, :] += node.y

        # shift arrow to node and plot
        ax.plot(node.x + rl[0, s:e], node.y + rl[1, s:e], 'k-')
        ax.add_patch(Polygon(np.transpose(rp), facecolor='k'))

    def plot_imposed_disp(self, ax, imposed_disp, max_disp):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        val = imposed_disp["val"] / max_disp
        node = imposed_disp["node"]
        dir = imposed_disp["dir"]

        small = self.get_small
        offset = 0.5 * small

        lf = abs(val) * 1.5 * small  # arrow length
        lh = 0.6 * small  # arrow head length
        wh = 0.6 * small  # arrow head width
        sp = 0.15 * small  # half spacing between double line
        lf = max(lf, lh * 1.5)

        (angle, num_el) = self.get_support_angle(node)
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        n = np.array([c, s])
        inward = (n[dir-1] == 0 or np.sign(n[dir-1]) == np.sign(val))

        to_rotate = (dir - 1) * 90 + (n[dir-1] >= 0) * 180
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
        rp[0, :] += node.x
        rp[1, :] += node.y

        ax.plot(rl[0, :] + node.x, rl[1, :] + node.y, 'k-')
        ax.add_patch(Polygon(np.transpose(rp),
                             facecolor='none', linewidth=1, edgecolor='k'))

    def plot_imposed_rot(self, ax, imposed_disp):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        node = imposed_disp["node"]
        val = imposed_disp["val"]

        small = self.get_small
        lh = 0.4 * small  # arrow head length
        wh = 0.4 * small  # arrow head width
        r1 = 1.0 * small
        r2 = 1.2 * small
        (angle, num_el) = self.get_support_angle(node)
        ths = np.arange(100, 261)

        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)
        rot_mat = np.array([[c, -s], [s, c]])

        # make arrow tail around (0,0)
        rr = (r1 + r2) / 2
        ll = np.array([rr * np.cos(ths * np.pi / 180),
                       rr * np.sin(ths * np.pi / 180)])
        l1 = np.array([r1 * np.cos(ths * np.pi / 180),
                       r1 * np.sin(ths * np.pi / 180)])
        l2 = np.array([r2 * np.cos(ths * np.pi / 180),
                       r2 * np.sin(ths * np.pi / 180)])

        # make arrow head at (0,0)
        pp = np.array([[-lh, -lh, 0], [-wh / 2, wh / 2, 0]])

        # rotate arrow head around (0,0)
        if val > 0:
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
        rp[0, :] += node.x
        rp[1, :] += node.y

        # shift arrow to node and plot
        ax.plot(node.x + rl1[0, :], node.y + rl1[1, :], 'k-')
        ax.plot(node.x + rl2[0, :], node.y + rl2[1, :], 'k-')
        ax.plot(node.x + np.append(rl1[0, ibase], rl2[0, ibase]),
                node.y + np.append(rl1[1, ibase], rl2[1, ibase]), 'k-')
        ax.add_patch(Polygon(np.transpose(rp),
                             facecolor='none', linewidth=1, edgecolor='k'))

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
