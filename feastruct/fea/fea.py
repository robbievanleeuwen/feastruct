import sys
import numpy as np
from matplotlib.patches import Polygon
from fea.exceptions import FEAInputError


class fea:
    """Parent class for a finite element analysis.

    Establishes a template for each different type of finite element analysis,
    e.g. frame, membrane etc. and provides a number of generic methods which
    are useful for all types of analyses.

    Attributes:
        nodes:          list of 'Node' items with which to initialise the
                        class.
        elements:       list of 'FiniteElement' items with which to initialise
                        the class.
        freedom_cases:  list of 'FreedomCase' itmes with which to initialise
                        the class.
        load_cases:     list of 'LoadCase' itmes with which to initialise
                        the class.
        analysis_cases: list of 'AnalysisCase' items with which to initialise
                        the class.
        non_linear:     boolean defining the type of analysis
                        (geometric & material).
        dofs:           structural degrees of freedom per node for current
                        analysis type.
    """

    def __init__(self, nodes, elements, freedom_cases, load_cases,
                 analysis_cases, non_linear, dofs):
        """Inits the fea class with..."""

        if nodes is None:
            self.nodes = []
        else:
            self.nodes = nodes

        if elements is None:
            self.elements = []
        else:
            self.elements = elements

        if freedom_cases is None:
            self.freedom_cases = []
        else:
            self.freedom_cases = freedom_cases

        if load_cases is None:
            self.load_cases = []
        else:
            self.load_cases = load_cases

        if analysis_cases is None:
            self.analysis_cases = []
        else:
            self.analysis_cases = analysis_cases

        self.non_linear = non_linear
        self.dofs = dofs

    def add_node(self, id, coord):
        """Adds a node to the fea class.

        Adds a Node object to the fea object. Checks for an exisiting id and
        ensures that a node is not already in the desired location.

        Args:
            id:     unique integer identifying the node
            coord:  cartesian coordinates of the node i.e. [x, y]

        Returns:
            void

        Raises:
            FEAInputError:  Raised if the id already exists or if a node
                            already exists at the desired location.
        """

        # TODO: check that node id and location does not already exist
        self.nodes.append(Node(id, coord))

        # raise exception if duplicate added

    def add_element(self, element):
        """Adds an element to the fea class.

        Adds a FiniteElement object to the fea object. Checks for an exisiting
        id.

        Args:
            element:    FiniteElement object to be added to the fea object.

        Returns:
            void

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check that element id does not already exist
        self.elements.append(element)

        # raise exception if duplicate added

    def add_freedom_case(self, id, items=None):
        """Adds a freedom case to the fea class.

        Adds a FreedomCase object to the fea object. Checks for an exisiting
        id.

        Args:
            id:             FreedomCase id.
            items:          List of items to initialise the FreedomCase with.

        Returns:
            fc:             FreedomCase object.

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check to see if id already exists and raise exception

        fc = FreedomCase(self, id, items)
        self.freedom_cases.append(fc)
        return fc

    def add_load_case(self, id, items=None):
        """Adds a load case to the fea class.

        Adds a LoadCase object to the fea object. Checks for an exisiting id.

        Args:
            id:             LoadCase id.
            items:          List of items to initialise the LoadCase with.

        Returns:
            lc:             LoadCase object.

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check to see if id already exists and raise exception

        lc = LoadCase(self, id, items)
        self.load_cases.append(lc)
        return lc

    def add_analysis_case(self, id, fc_id, lc_id):
        """Adds an analysis case to the fea class.

        Adds an AnalysisCase object to the fea object. Checks for an exisiting
        id.

        Args:
            id:             AnalysisCase id.
            fc_id:          FreedomCase id to add to the analysis case.
            lc_id:          LoadCase id to add to the analysis case.

        Returns:
            void

        Raises:
            FEAInputError:  Raised if the id already exists.
        """

        # TODO: check to see if id already exists and raise exception

        self.analysis_cases.append(AnalysisCase(self, id, fc_id, lc_id))

    def find_node(self, node_id):
        """Finds the reference to the node object given its node_id.
        """

        # return the Node object whose id matches the given node_id
        for node in self.nodes:
            if node.id == node_id:
                return node
        else:
            raise FEAInputError("Cannot find node_id: {}".format(node_id))

    def get_node_lims(self):
        """asldkjasld
        """

        xmin = self.nodes[0].x
        xmax = self.nodes[0].x
        ymin = self.nodes[0].y
        ymax = self.nodes[0].y

        for node in self.nodes:
            xmin = min(xmin, node.x)
            xmax = max(xmax, node.x)
            ymin = min(ymin, node.y)
            ymax = max(ymax, node.y)

        return (xmin, xmax, ymin, ymax)

    def find_analysis_case(self, case_id):
        """
        """

        # return the AnalysisCase object whose id matches the given case_id
        for analysis_case in self.analysis_cases:
            if analysis_case.id == case_id:
                return analysis_case
        else:
            raise FEAInputError("Cannot find AnalysisCase id: {}".format(
                case_id))


class FiniteElement:
    """
    """

    def __init__(self, analysis, id, node_ids):
        """Inits the FiniteElement class with an analysis object, element id
        and corresponding node ids.

        Args:
            analysis: reference to the analysis object.
            id: an integer representing a unique element id.
            node_ids: a list containing the node ids defining the element.

        Raises:
            TypeError: TODO
            ValueError: Raised if a negative id is provided. TODO
        """

        self.analysis = analysis
        self.id = id
        self.nodes = []
        self.f_int = []

        # add references to the node objects
        self.get_nodes(node_ids)

    def get_nodes(self, node_ids):
        """
        """

        for node_id in node_ids:
            try:
                self.nodes.append(self.analysis.find_node(node_id))
            except FEAInputError as error:
                print(
                    "Error in FiniteElement id: {}. {}".format(self.id, error))
                sys.exit(1)

    def get_node_coords(self):
        """
        """

        coord_list = []

        for node in self.nodes:
            coord_list.append(node.coords)

        return np.array(coord_list)

    def get_nodal_displacements(self, case_id):
        """
        """

        disp_list = []  # allocate list of displacements

        for node in self.nodes:
            disp_list.append(node.get_displacement(case_id))

        return np.array(disp_list)

    def get_dofs(self):
        """
        """

        dof_list = np.array([], dtype=np.int32)

        for node in self.nodes:
            dof_list = np.append(dof_list, node.dofs)

        return dof_list

    def get_stiff_matrix(self):
        """
        """

        pass

    def save_fint(self, f_int, case_id):
        """
        """

        # save internal force vector in global coordinate system
        self.f_int.append({"case_id": case_id, "f_int": f_int})


class Node:
    """
    """

    def __init__(self, id, coord):
        """Inits the Node class with an id and coordinate.

        Args:
            id: an integer representing a unique node id.
            coord: a list consisting of the x, y coordinates of
            the node.

        Raises:
            TypeError: TODO
            ValueError: Raised if a negative id is provided. TODO
        """

        self.id = id
        self.coord = coord
        self.dofs = []
        self.u = []
        self.fixity = [0, 0, 0]  # for post processing only
        # TODO: check value types

    @property
    def x(self):
        return self.coord[0]

    @property
    def y(self):
        return self.coord[1]

    @property
    def coords(self):
        return [self.x, self.y]

    def get_displacement(self, case_id):
        """
        """

        # get dictionary displacement entry for given case_id
        disp = next(d for d in self.u if d["case_id"] == case_id)

        if disp is not None:
            return disp["u"]
        else:
            raise FEAInputError("""Cannot find an analysis result for
            case_id: {} at node_id: {}""".format(case_id, self.id))


class Case:
    """
    """

    def __init__(self, analysis, id, items):
        """"""

        self.analysis = analysis
        self.id = id

        if items is None:
            self.items = []
        else:
            self.items = items

    def add_item(self, item):
        """"""

        self.items.append(item)


class FreedomCase(Case):
    """
    """

    def __init__(self, analysis, id, items):
        super().__init__(analysis, id, items)

    def add_nodal_support(self, node_id, val, dir):
        """
        """

        # TODO: check that the support does not already exist
        # raise exception if duplicate added

        # add a dictionary entry to the supports list
        self.add_item(NodalSupport(self.analysis, node_id, val, dir))


class LoadCase(Case):
    """
    """

    def __init__(self, analysis, id, items):
        super().__init__(analysis, id, items)

    def add_nodal_load(self, node_id, val, dir):
        """
        """

        # TODO: check that the load does not already exist
        # raise exception if duplicate added

        # add a dictionary entry to the supports list
        self.add_item(NodalLoad(self.analysis, node_id, val, dir))


class AnalysisCase:
    """
    """

    def __init__(self, analysis, id, fc_id, lc_id):
        """
        """

        self.analysis = analysis
        self.id = id

        try:
            # find freedom case
            for fc in self.analysis.freedom_cases:
                if fc.id == fc_id:
                    self.freedom_case = fc
                    break
            else:
                raise FEAInputError("Cannot find FreedomCase id: {}".format(
                    fc_id))

            # find load case
            for lc in self.analysis.load_cases:
                if lc.id == lc_id:
                    self.load_case = lc
                    break
            else:
                raise FEAInputError("Cannot find LoadCase id: {}".format(
                    lc_id))
        except FEAInputError as error:
            print(error)
            sys.exit(1)


class BoundaryCondition:
    """asdlkjaskld
    """

    def __init__(self, analysis, node_id, val, dir):
        # TODO: check value types e.g. node_id and dir are ints, check dir is
        # within dofs limits

        # find the node object corresponding to node_id
        try:
            node = analysis.find_node(node_id)
        except IndexError as error:
            print(error)
            sys.exit(1)

        self.node = node
        self.val = val
        self.dir = dir


class NodalSupport(BoundaryCondition):
    """asldkjasdl
    """

    def __init__(self, analysis, node_id, val, dir):
        super().__init__(analysis, node_id, val, dir)

        self.reaction = []

    def apply_support(self, K, f_ext):
        """sadsad
        """

        # modify stiffness matrix and f_ext
        K[self.node.dofs[self.dir-1], :] = 0
        K[self.node.dofs[self.dir-1], self.node.dofs[self.dir-1]] = 1
        f_ext[self.node.dofs[self.dir-1]] = self.val

    def plot_support(self, ax, max_disp, small, get_support_angle, case_id,
                     deformed, def_scale):
        """asdasdas
        """

        if self.node.fixity is not [1, 1, 0]:
            (angle, num_el) = get_support_angle(self.node)

        if self.node.fixity == [1, 0, 0]:
            # ploy a y-roller
            angle = round(angle / 180) * 180
            self.plot_xysupport(ax, angle, True, num_el == 1, small, case_id,
                                deformed, def_scale)

        elif self.node.fixity == [0, 1, 0]:
            # plot an x-roller
            if np.mod(angle + 1, 180) < 2:  # prefer support below
                angle = 90
            else:
                angle = round((angle + 90) / 180) * 180 - 90

            self.plot_xysupport(ax, angle, True, num_el == 1, small, case_id,
                                deformed, def_scale)

        elif self.node.fixity == [1, 1, 0]:
            # plot a hinge
            (angle, num_el) = get_support_angle(self.node, 2)
            self.plot_xysupport(ax, angle, False, num_el == 1, small, case_id,
                                deformed, def_scale)

        elif self.node.fixity == [0, 0, 1]:
            ax.plot(self.node.x, self.node.y, 'kx', markersize=8)

        else:
            # plot a support with moment fixity
            if self.node.fixity == [1, 1, 1]:
                # plot a fixed support
                s = np.sin(angle * np.pi / 180)
                c = np.cos(angle * np.pi / 180)
                rot_mat = np.array([[c, -s], [s, c]])
                line = np.array([[0, 0], [-1, 1]]) * small
                rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * small
                ec = 'none'

            elif self.node.fixity == [1, 0, 1]:
                # plot y-roller block
                angle = round(angle / 180) * 180
                s = np.sin(angle * np.pi / 180)
                c = np.cos(angle * np.pi / 180)
                rot_mat = np.array([[c, -s], [s, c]])
                line = np.array([[-0.85, -0.85], [-1, 1]]) * small
                rect = np.array([[-0.6, -0.6, 0, 0], [-1, 1, 1, -1]]) * small
                ec = 'k'

            elif self.node.fixity == [0, 1, 1]:
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
                u = self.node.get_displacement(case_id)

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
            ax.add_patch(Polygon(np.transpose(rot_rect),
                                 facecolor=(0.7, 0.7, 0.7), edgecolor=ec))

    def plot_xysupport(self, ax, angle, roller, hinge, small, case_id,
                       deformed, def_scale):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        # determine coordinates of node
        if deformed:
            # get displacement of node for current analysis case
            u = self.node.get_displacement(case_id)

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
            ax.add_patch(Polygon(np.transpose(rot_rect),
                                 facecolor=(0.7, 0.7, 0.7)))

        ax.plot(rot_triangle[0, :] + x, rot_triangle[1, :] + y, 'k-',
                linewidth=1)

        if hinge:
            ax.plot(x, y, 'ko', markerfacecolor='w', linewidth=1, markersize=4)

    def plot_imposed_disp(self, ax, max_disp, small, get_support_angle,
                          case_id, deformed, def_scale):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
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
        inward = (n[self.dir-1] == 0 or np.sign(n[self.dir-1]) == np.sign(val))

        to_rotate = (self.dir - 1) * 90 + (n[self.dir-1] >= 0) * 180
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
            u = self.node.get_displacement(case_id)

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
        ax.add_patch(Polygon(np.transpose(rp),
                             facecolor='none', linewidth=1, edgecolor='k'))

    def plot_imposed_rot(self, ax, small, get_support_angle, case_id, deformed,
                         def_scale):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
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
        ll = np.array([rr * np.cos(ths * np.pi / 180),
                       rr * np.sin(ths * np.pi / 180)])
        l1 = np.array([r1 * np.cos(ths * np.pi / 180),
                       r1 * np.sin(ths * np.pi / 180)])
        l2 = np.array([r2 * np.cos(ths * np.pi / 180),
                       r2 * np.sin(ths * np.pi / 180)])

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
            u = self.node.get_displacement(case_id)

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
        ax.plot(np.append(rl1[0, ibase], rl2[0, ibase]),
                np.append(rl1[1, ibase], rl2[1, ibase]), 'k-')
        ax.add_patch(Polygon(np.transpose(rp),
                             facecolor='none', linewidth=1, edgecolor='k'))


class NodalLoad(BoundaryCondition):
    """asldkjasdl
    """

    def __init__(self, analysis, node_id, val, dir):
        super().__init__(analysis, node_id, val, dir)

    def apply_load(self, f_ext):
        """alskdjaskld
        """

        # add load to f_ext, selecting the correct dof from dofs
        f_ext[self.node.dofs[self.dir-1]] = self.val

    def plot_load(self, ax, max_force, small, get_support_angle, case_id,
                  deformed, def_scale):
        """aslkdjsak

        N.B. this method is adopted from the MATLAB code by F.P. van der Meer:
        plotGeom.m.
        """

        # determine coordinates of node
        if deformed:
            # get displacement of node for current analysis case
            u = self.node.get_displacement(case_id)

            x = self.node.x + u[0] * def_scale
            y = self.node.y + u[1] * def_scale
        else:
            x = self.node.x
            y = self.node.y

        val = self.val / max_force

        offset = 0.5 * small
        (angle, num_el) = get_support_angle(self.node)
        s = np.sin(angle * np.pi / 180)
        c = np.cos(angle * np.pi / 180)

        # plot nodal force
        if self.dir not in (3, 6):
            lf = abs(val) * 1.5 * small  # arrow length
            lh = 0.6 * small  # arrow head length
            wh = 0.6 * small  # arrow head width
            lf = max(lf, lh * 1.5)

            n = np.array([c, s])
            inward = (n[self.dir-1] == 0 or
                      np.sign(n[self.dir-1]) == np.sign(val))

            to_rotate = (self.dir - 1) * 90 + (n[self.dir-1] > 0) * 180
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
            ll = np.array([rr * np.cos(ths * np.pi / 180),
                           rr * np.sin(ths * np.pi / 180)])

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
